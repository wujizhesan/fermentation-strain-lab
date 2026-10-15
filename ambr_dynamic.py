import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

import numpy as np
import openpyxl


ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT/'research/ambr_dynamic_protocol.json'
CHANNELS = {'CER':7,'OUR':14,'CO2_pct':8,'O2_pct':13,'DO':9,'pH':15,
            'feed_rate':10,'feed_volume':11,'temperature':19,'air_flow':3,'stir_speed':18}
RECIPE = ['time_h','time_squared','strain_D161646','pH_setpoint','YE_nominal_g_L','glucose_feed_nominal_g_L']
HISTORY_CHANNELS = ['CER','OUR','CO2_pct','DO']
TARGETS = {'CDW_g_L':32,'mucic_acid_corrected_g_L':31}
STRATEGIES = ['training_mean','time_recipe_ridge','instantaneous_ridge','history_ridge']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def save_csv(path, rows):
    with Path(path).open('w',encoding='utf-8',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def validate_source(data, spec):
    if len(data)!=spec['bytes'] or hashlib.sha256(data).hexdigest()!=spec['sha256'] or hashlib.md5(data).hexdigest()!=spec['published_md5']:
        raise ValueError('Published source checksum mismatch')


def metadata(name):
    match = re.fullmatch(r'(D-161646|T2) (pH|YE|glucose) (\d+(?:\.\d+)?) \(([12])\)',name)
    if not match:
        raise ValueError('Unknown culture identity: '+name)
    strain,kind,value,round_id = match.groups()
    value = float(value)
    return {'batch_id':name,'round_id':int(round_id),'condition':name[:-4],
            'strain':'D-161646' if strain=='D-161646' else 'D-221704',
            'strain_D161646':float(strain=='D-161646'),
            'pH_setpoint':value if kind=='pH' else 4.,
            'YE_nominal_g_L':value if kind=='YE' else (1. if kind=='glucose' else 2.),
            'glucose_feed_nominal_g_L':value if kind=='glucose' else 30.}


def read_source(path):
    book = openpyxl.load_workbook(path,read_only=True,data_only=True)
    description = [str(r[0]) for r in book['Description'].values if r[0] is not None]
    feed = {r[0]:{'galacturonate_g_L':float(r[1]),'glucose_g_L':float(r[2])} for r in list(book['Feed composition'].values)[4:] if r[0]}
    cultures = []
    for sheet in book:
        if sheet.title in ['Description','Feed composition']:
            continue
        meta = metadata(sheet.title)
        cells = list(sheet.values)
        headers = cells[0]
        if headers[0]!='Time' or headers[24]!='Time (h)' or headers[31]!='Mucic acid (g/L)' or headers[32]!='CDW (g/L)':
            raise ValueError('Source schema changed')
        for key,index in CHANNELS.items():
            expected = {'CO2_pct':'CO2%','O2_pct':'O2%','feed_rate':'Feed#1 flow rate','feed_volume':'Feed#1 volume pumped',
                        'temperature':'Temperature','air_flow':'Air flow','stir_speed':'Stir speed'}.get(key,key)
            if headers[index].split('-',1)[1].strip()!=expected:
                raise ValueError('Online channel changed: '+key)
        online,offline = [],[]
        for source_row,values in enumerate(cells[1:],2):
            if isinstance(values[0],(float,int)):
                vector = [values[0]]+[values[j] for j in CHANNELS.values()]
                if not all(isinstance(v,(int,float)) and np.isfinite(v) for v in vector):
                    raise ValueError('Missing selected online signal')
                online.append([source_row]+[float(v) for v in vector])
            if isinstance(values[24],(float,int)):
                r = {**meta,'row_id':sheet.title+'!'+str(source_row),'source_sheet':sheet.title,'source_row':source_row,'time_h':float(values[24])}
                for target,index in TARGETS.items():
                    value = values[index]
                    if value is not None and (not isinstance(value,(float,int)) or not np.isfinite(value) or value<0):
                        raise ValueError('Invalid measured target')
                    r[target] = float(value) if value is not None else None
                offline.append(r)
        online = np.asarray(online,dtype=float)
        if not len(online) or np.any(np.diff(online[:,1])<=0):
            raise ValueError('Online clock must be strictly increasing')
        if len({r['time_h'] for r in offline})!=len(offline) or any(a['time_h']>b['time_h'] for a,b in zip(offline,offline[1:])):
            raise ValueError('Offline sampling identities are ambiguous')
        cultures.append({'metadata':meta,'online':online,'offline':offline,'feed_composition':feed[sheet.title]})
    book.close()
    return cultures,description


def feature_names(strategy):
    if strategy=='training_mean':
        return []
    if strategy not in STRATEGIES:
        raise ValueError('Unknown strategy')
    result = list(RECIPE)
    if strategy!='time_recipe_ridge':
        result += ['last_'+c for c in CHANNELS]
    if strategy=='history_ridge':
        result += [stat+'_'+c+('_'+str(int(window))+'h' if stat=='mean' else '')
                   for c in HISTORY_CHANNELS for stat,window in [('mean',1.),('mean',6.),('integral',0.)]]
    return result


def snapshot(online, row, max_age_h=.5):
    query = float(row['time_h'])
    if not np.isfinite(query) or query<0 or max_age_h<0:
        raise ValueError('Invalid sampling time or age threshold')
    if len(online) and np.any(np.diff(online[:,1])<=0):
        raise ValueError('Invalid online ordering')
    past = online[online[:,1]<=query]
    if not len(past):
        return None,'no_past_online_observation'
    age = query-past[-1,1]
    if age>max_age_h:
        return None,'online_observation_too_old'
    features = {'time_h':query,'time_squared':query**2,**{k:float(row[k]) for k in RECIPE[2:]}}
    channel_columns = {c:i+2 for i,c in enumerate(CHANNELS)}
    for channel,index in channel_columns.items():
        features['last_'+channel] = float(past[-1,index])
    for channel in HISTORY_CHANNELS:
        values = past[:,channel_columns[channel]]
        for hours in [1.,6.]:
            window = past[:,1]>=query-hours
            features['mean_'+channel+'_'+str(int(hours))+'h'] = float(values[window].mean())
        features['integral_'+channel] = float(np.trapezoid(values,past[:,1]))
    if not np.all(np.isfinite(list(features.values()))):
        raise ValueError('Nonfinite causal feature')
    return {**features,'latest_online_source_row':int(past[-1,0]),'latest_online_time_h':float(past[-1,1]),
            'online_age_h':float(age),'prefix_points':len(past)},''


def partitions(rows, train_round, target):
    if train_round not in [1,2] or target not in TARGETS:
        raise ValueError('Unsupported split or target')
    eligible = [r for r in rows if r['feature_eligible'] and r[target] is not None]
    train = [r for r in eligible if r['round_id']==train_round]
    test = [r for r in eligible if r['round_id']!=train_round]
    if not train or not test or {r['batch_id'] for r in train}&{r['batch_id'] for r in test}:
        raise ValueError('Invalid independent-round split')
    return train,test


def fit(rows, target, strategy, alpha):
    names = feature_names(strategy)
    y = np.asarray([r[target] for r in rows])
    if not len(rows) or np.any(~np.isfinite(y)) or not np.isfinite(alpha) or alpha<=0:
        raise ValueError('Invalid training data or ridge penalty')
    model = {'strategy':strategy,'target':target,'features':names,'alpha':alpha,'intercept':float(y.mean()),
             'train_row_ids':[r['row_id'] for r in rows]}
    if names:
        x = np.asarray([[r[n] for n in names] for r in rows])
        center,scale = x.mean(axis=0),x.std(axis=0)
        scale[scale==0] = 1.
        z = (x-center)/scale
        coefficients = np.linalg.solve(z.T@z+alpha*np.eye(z.shape[1]),z.T@(y-y.mean()))
        model.update({'center':center.tolist(),'scale':scale.tolist(),'coefficients':coefficients.tolist()})
    return model


def predict(model, rows):
    names = model['features']
    if not names:
        return np.full(len(rows),model['intercept'])
    x = np.asarray([[r[n] for n in names] for r in rows])
    return model['intercept']+((x-model['center'])/model['scale'])@np.asarray(model['coefficients'])


def metrics(rows):
    error = np.asarray([r['prediction']-r['truth'] for r in rows])
    return {'points':len(rows),'RMSE':float(np.sqrt(np.mean(error**2))),'MAE':float(np.mean(np.abs(error)))}


def run(output, download):
    protocol = json.loads(PROTOCOL.read_text(encoding='utf-8'))
    if not protocol.get('locked_at_utc'):
        raise ValueError('Freeze protocol before fitting')
    for name,expected in protocol['pinned_inputs'].items():
        if digest(ROOT/name)!=expected:
            raise ValueError('Frozen input changed: '+name)
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use an empty output directory')
    spec = protocol['source_file']
    path = ROOT/spec['path']
    if not path.exists():
        if not download:
            raise FileNotFoundError(path)
        with urlopen(Request(spec['url'],headers={'User-Agent':'fermentation-research/1.0'}),timeout=60) as response:
            data = response.read(spec['bytes']+1)
        validate_source(data,spec)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(data)
    validate_source(path.read_bytes(),spec)
    cultures,description = read_source(path)
    rows,inventory,channel_audit = [],[],[]
    for culture in cultures:
        online,offline = culture['online'],culture['offline']
        inventory.append({**culture['metadata'],'online_points':len(online),'offline_points':len(offline),
                          'first_online_time_h':float(online[0,1]),'last_online_time_h':float(online[-1,1]),
                          'median_interval_h':float(np.median(np.diff(online[:,1]))),**culture['feed_composition']})
        channel_audit.extend({'batch_id':culture['metadata']['batch_id'],'channel':c,'points':len(online),
                              'minimum':float(online[:,i+2].min()),'maximum':float(online[:,i+2].max()),
                              'zero_points':int(np.count_nonzero(online[:,i+2]==0))} for i,c in enumerate(CHANNELS))
        for r in offline:
            features,reason = snapshot(online,r,protocol['maximum_online_age_h'])
            rows.append({**r,'feature_eligible':features is not None,'feature_issue':reason,
                         **(features or {n:None for n in feature_names('history_ridge')[len(RECIPE):]+['time_squared','latest_online_source_row','latest_online_time_h','online_age_h','prefix_points']})})
    output.mkdir(parents=True,exist_ok=True)
    (output/'protocol.json').write_bytes(PROTOCOL.read_bytes())
    save_csv(output/'observations.csv',rows)
    save_csv(output/'batch_inventory.csv',inventory)
    save_csv(output/'online_channel_audit.csv',channel_audit)
    models,predictions,folds,summary = [],[],[],[]
    for train_round in protocol['training_rounds']:
        direction = str(train_round)+'_to_'+str(3-train_round)
        for target in TARGETS:
            train,test = partitions(rows,train_round,target)
            for strategy in STRATEGIES:
                model = fit(train,target,strategy,protocol['ridge_alpha'])
                model.update({'model_id':len(models),'direction':direction,'train_round':train_round,'test_round':3-train_round})
                models.append(model)
                raw = predict(model,test)
                scored = [{ 'row_id':r['row_id'],'batch_id':r['batch_id'],'round_id':r['round_id'],'condition':r['condition'],
                            'time_h':r['time_h'],'direction':direction,'target':target,'strategy':strategy,'model_id':model['model_id'],
                            'truth':r[target],'raw_prediction':float(p),'prediction':float(max(0.,p))} for r,p in zip(test,raw)]
                predictions.extend(scored)
                groups = defaultdict(list)
                for r in scored:
                    groups[r['batch_id']].append(r)
                for batch,points in groups.items():
                    folds.append({'direction':direction,'target':target,'strategy':strategy,'batch_id':batch,**metrics(points)})
                ff = [f for f in folds if f['direction']==direction and f['target']==target and f['strategy']==strategy]
                summary.append({'direction':direction,'target':target,'strategy':strategy,'batches':len(ff),
                                'batch_macro_RMSE':float(np.mean([r['RMSE'] for r in ff])),
                                'batch_macro_MAE':float(np.mean([r['MAE'] for r in ff])),
                                'negative_raw_predictions':int(np.sum(raw<0)),**{'pooled_'+k:v for k,v in metrics(scored).items()}})
    pairs = []
    by_fold = {(r['direction'],r['target'],r['strategy'],r['batch_id']):r for r in folds}
    for r in folds:
        if r['strategy']!='history_ridge':
            continue
        baseline = by_fold[(r['direction'],r['target'],'instantaneous_ridge',r['batch_id'])]
        pairs.append({'direction':r['direction'],'target':r['target'],'batch_id':r['batch_id'],
                      'history_RMSE':r['RMSE'],'instantaneous_RMSE':baseline['RMSE'],'delta_RMSE':r['RMSE']-baseline['RMSE']})
    save_json(output/'models.json',models)
    save_csv(output/'predictions.csv',predictions)
    save_csv(output/'fold_metrics.csv',folds)
    save_csv(output/'summary.csv',summary)
    save_csv(output/'paired_batches.csv',pairs)
    artifacts = ['protocol.json','observations.csv','batch_inventory.csv','online_channel_audit.csv','models.json','predictions.csv','fold_metrics.csv','summary.csv','paired_batches.csv']
    primary = [r for r in pairs if r['direction']=='1_to_2' and r['target']=='CDW_g_L']
    report = {'created_at_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':spec['sha256'],
              'protocol_sha256':digest(PROTOCOL),'implementation_sha256':digest(__file__),
              'cultures':len(cultures),'independent_rounds':len({r['round_id'] for r in rows}),
              'online_rows':sum(len(c['online']) for c in cultures),'offline_rows':len(rows),
              'feature_eligible_rows':sum(r['feature_eligible'] for r in rows),'feature_excluded_rows':[{'row_id':r['row_id'],'reason':r['feature_issue']} for r in rows if not r['feature_eligible']],
              'eligible_points_by_target':{t:sum(r['feature_eligible'] and r[t] is not None for r in rows) for t in TARGETS},
              'model_fits':len(models),'prediction_rows':len(predictions),'fold_scores':len(folds),'summaries':summary,
              'primary_comparison':{'direction':'1_to_2','target':'CDW_g_L','candidate':'history_ridge','reference':'instantaneous_ridge',
                                    'batches':len(primary),'mean_batch_RMSE_delta':float(np.mean([r['delta_RMSE'] for r in primary])),
                                    'batches_improved':sum(r['delta_RMSE']<0 for r in primary)},
              'workbook_description':description,'scope_limits':protocol['scope_limits'],
              'artifact_sha256':{n:digest(output/n) for n in artifacts}}
    save_json(output/'report.json',report)
    print(json.dumps({k:report[k] for k in ['cultures','online_rows','offline_rows','eligible_points_by_target','model_fits','prediction_rows','primary_comparison']},indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'research/results/ambr_dynamic')
    parser.add_argument('--download',action='store_true')
    args = parser.parse_args()
    run(args.output,args.download)
