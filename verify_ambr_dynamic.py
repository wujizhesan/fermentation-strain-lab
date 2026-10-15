import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

import numpy as np


ROOT = Path(__file__).resolve().parent
NAMESPACE = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
GASES = {'CER':'H','OUR':'O','CO2_pct':'I','O2_pct':'N','DO':'J','pH':'P',
         'feed_rate':'K','feed_volume':'L','temperature':'T','air_flow':'D','stir_speed':'S'}
BASE = ['time_h','time_squared','strain_D161646','pH_setpoint','YE_nominal_g_L','glucose_feed_nominal_g_L']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def read_csv(path):
    with Path(path).open(encoding='utf-8',newline='') as stream:
        return list(csv.DictReader(stream))


def raw_workbook(path):
    output = {}
    with zipfile.ZipFile(path) as archive:
        shared = ET.fromstring(archive.read('xl/sharedStrings.xml')) if 'xl/sharedStrings.xml' in archive.namelist() else None
        strings = [''.join(t.itertext()) for t in shared.findall('s:si',NAMESPACE)] if shared is not None else []
        rels = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        targets = {r.get('Id'):r.get('Target') for r in rels}
        doc = ET.fromstring(archive.read('xl/workbook.xml'))
        for sheet in doc.findall('s:sheets/s:sheet',NAMESPACE):
            title = sheet.get('name')
            ref = sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
            target = targets[ref]
            member = target.lstrip('/') if target.startswith('/') else 'xl/'+target
            tree = ET.fromstring(archive.read(member))
            rows = []
            for row in tree.findall('s:sheetData/s:row',NAMESPACE):
                cells = {}
                for cell in row.findall('s:c',NAMESPACE):
                    column = re.match(r'[A-Z]+',cell.get('r')).group()
                    value = cell.find('s:v',NAMESPACE)
                    kind = cell.get('t')
                    if kind=='inlineStr':
                        cells[column] = ''.join(cell.find('s:is',NAMESPACE).itertext())
                    elif value is not None and value.text:
                        cells[column] = strings[int(value.text)] if kind=='s' else float(value.text)
                rows.append((int(row.get('r')),cells))
            output[title] = rows
    return output


def verify(output):
    protocol = read_json(output/'protocol.json')
    report = read_json(output/'report.json')
    assert protocol==read_json(ROOT/'research/ambr_dynamic_protocol.json')
    assert sha(output/'protocol.json')==report['protocol_sha256']
    assert sha(ROOT/'ambr_dynamic.py')==report['implementation_sha256']
    for name,expected in protocol['pinned_inputs'].items():
        assert sha(ROOT/name)==expected,name
    for name,expected in report['artifact_sha256'].items():
        assert sha(output/name)==expected,name
    spec = protocol['source_file']
    source = ROOT/spec['path']
    assert source.stat().st_size==spec['bytes']
    assert sha(source)==spec['sha256']==report['source_sha256']
    assert hashlib.md5(source.read_bytes()).hexdigest()==spec['published_md5']
    raw = raw_workbook(source)
    prepared = {r['row_id']:r for r in read_csv(output/'observations.csv')}
    independent,inventory,raw_channels = {},{},{}
    feed = {r['A']:(r['B'],r['C']) for i,r in raw['Feed composition'] if i>=5 and r.get('A')}
    feature_error = 0.
    online_points,offline_points,excluded = 0,0,0
    feature_lists = {'training_mean':[], 'time_recipe_ridge':BASE,
                     'instantaneous_ridge':BASE+['last_'+c for c in GASES]}
    feature_lists['history_ridge'] = feature_lists['instantaneous_ridge']+[x for c in ['CER','OUR','CO2_pct','DO'] for x in ['mean_'+c+'_1h','mean_'+c+'_6h','integral_'+c]]
    for name,rows in raw.items():
        if name in ['Description','Feed composition']:
            continue
        strain,kind,number,replicate = re.fullmatch(r'(D-161646|T2) (pH|YE|glucose) (\d+(?:\.\d+)?) \(([12])\)',name).groups()
        round_id = int(replicate)
        number = float(number)
        online = [(i,r) for i,r in rows[1:] if isinstance(r.get('A'),(float,int))]
        samples = [(i,r) for i,r in rows[1:] if isinstance(r.get('Y'),(float,int))]
        online_points += len(online)
        offline_points += len(samples)
        for channel,column in GASES.items():
            signal = [r[column] for _,r in online]
            raw_channels[(name,channel)] = (len(signal),min(signal),max(signal),sum(v==0 for v in signal))
        times = [r['A'] for _,r in online]
        assert all(a<b for a,b in zip(times,times[1:]))
        inventory[name] = (len(online),len(samples),times[0],times[-1])
        for source_row,values in samples:
            row_id = name+'!'+str(source_row)
            observed = prepared[row_id]
            time = values['Y']
            assert observed['source_sheet']==observed['batch_id']==name and observed['condition']==name[:-4]
            assert int(observed['source_row'])==source_row and int(observed['round_id'])==round_id
            assert observed['strain']==('D-161646' if strain=='D-161646' else 'D-221704')
            assert float(observed['time_h'])==time
            features = {'time_h':time,'time_squared':time*time,'strain_D161646':float(strain=='D-161646'),
                        'pH_setpoint':number if kind=='pH' else 4.,'YE_nominal_g_L':number if kind=='YE' else (1. if kind=='glucose' else 2.),
                        'glucose_feed_nominal_g_L':number if kind=='glucose' else 30.}
            for target,column in [('CDW_g_L','AG'),('mucic_acid_corrected_g_L','AF')]:
                actual = values.get(column)
                assert (observed[target]=='')==(actual is None)
                if actual is not None:
                    assert float(observed[target])==actual and actual>=0
            prefix = [(i,r) for i,r in online if r['A']<=time]
            issue = 'no_past_online_observation' if not prefix else 'online_observation_too_old' if time-prefix[-1][1]['A']>protocol['maximum_online_age_h'] else ''
            assert (observed['feature_eligible']=='True')==(issue=='')
            assert observed['feature_issue']==issue
            if issue:
                excluded += 1
                continue
            latest_row,latest = prefix[-1]
            assert int(observed['latest_online_source_row'])==latest_row
            assert float(observed['latest_online_time_h'])==latest['A']<=time
            assert float(observed['online_age_h'])==time-latest['A']
            assert int(observed['prefix_points'])==len(prefix)
            for channel,column in GASES.items():
                features['last_'+channel] = latest[column]
            for channel in ['CER','OUR','CO2_pct','DO']:
                column = GASES[channel]
                for hours in [1,6]:
                    signal = [r[column] for _,r in prefix if r['A']>=time-hours]
                    features['mean_'+channel+'_'+str(hours)+'h'] = sum(signal)/len(signal)
                features['integral_'+channel] = sum((a[column]+b[column])*(b['A']-a['A'])/2 for (_,a),(_,b) in zip(prefix,prefix[1:]))
            for feature,value in features.items():
                actual = float(observed[feature])
                feature_error = max(feature_error,abs(actual-value))
                assert np.isclose(actual,value,rtol=1e-10,atol=1e-9),(row_id,feature,actual,value)
            independent[row_id] = {'features':features,'round_id':round_id,'batch_id':name,
                                   'CDW_g_L':values.get('AG'),'mucic_acid_corrected_g_L':values.get('AF')}
    assert len(prepared)==offline_points==report['offline_rows']==230
    assert online_points==report['online_rows']==74849 and len(inventory)==report['cultures']==20
    assert excluded==len(report['feature_excluded_rows'])
    assert {r['row_id']:r['reason'] for r in report['feature_excluded_rows']}=={i:r['feature_issue'] for i,r in prepared.items() if r['feature_eligible']=='False'}
    assert report['feature_eligible_rows']==len(independent)
    assert report['independent_rounds']==len({r['round_id'] for r in independent.values()})==2
    assert report['eligible_points_by_target']=={t:sum(r[t] is not None for r in independent.values()) for t in ['CDW_g_L','mucic_acid_corrected_g_L']}
    assert report['workbook_description']==[r['A'] for _,r in raw['Description'] if 'A' in r]
    for r in read_csv(output/'online_channel_audit.csv'):
        assert (int(r['points']),float(r['minimum']),float(r['maximum']),int(r['zero_points']))==raw_channels[(r['batch_id'],r['channel'])]
    for r in read_csv(output/'batch_inventory.csv'):
        original = inventory[r['batch_id']]
        assert (int(r['online_points']),int(r['offline_points']),float(r['first_online_time_h']),float(r['last_online_time_h']))==original
        assert (float(r['galacturonate_g_L']),float(r['glucose_g_L']))==feed[r['batch_id']]
    models = read_json(output/'models.json')
    assert {(m['train_round'],m['target'],m['strategy']) for m in models}=={(replicate,t,s) for replicate in [1,2] for t in ['CDW_g_L','mucic_acid_corrected_g_L'] for s in feature_lists}
    predictions = read_csv(output/'predictions.csv')
    model_error,prediction_error = 0.,0.
    checked = 0
    for model in models:
        target,train_round,names = model['target'],model['train_round'],feature_lists[model['strategy']]
        assert model['features']==names and model['alpha']==protocol['ridge_alpha']==10.
        train = [(i,r) for i,r in independent.items() if r['round_id']==train_round and r[target] is not None]
        test = [(i,r) for i,r in independent.items() if r['round_id']==3-train_round and r[target] is not None]
        assert set(model['train_row_ids'])=={i for i,_ in train}
        assert not {r['batch_id'] for _,r in train}&{r['batch_id'] for _,r in test}
        y = np.array([r[target] for _,r in train])
        intercept = y.mean()
        assert np.isclose(intercept,model['intercept'],rtol=1e-12,atol=1e-10)
        if names:
            x = np.array([[r['features'][n] for n in names] for _,r in train])
            center,scale = np.mean(x,axis=0),np.std(x,axis=0)
            scale[scale==0] = 1.
            z = (x-center)/scale
            expanded = np.vstack([z,np.sqrt(model['alpha'])*np.eye(len(names))])
            coef = np.linalg.lstsq(expanded,np.concatenate([y-intercept,np.zeros(len(names))]),rcond=None)[0]
            for a,b in [(center,model['center']),(scale,model['scale']),(coef,model['coefficients'])]:
                model_error = max(model_error,float(np.max(np.abs(a-b))))
                assert np.allclose(a,b,rtol=1e-8,atol=1e-8)
            pp = intercept+(np.array([[r['features'][n] for n in names] for _,r in test])-center)/scale@coef
        else:
            pp = np.full(len(test),intercept)
        by_id = {r['row_id']:r for r in predictions if int(r['model_id'])==model['model_id']}
        assert set(by_id)=={i for i,_ in test}
        for (row_id,r),prediction in zip(test,pp):
            saved = by_id[row_id]
            assert saved['batch_id']==r['batch_id'] and saved['target']==target and saved['strategy']==model['strategy']
            assert saved['direction']==model['direction']==str(train_round)+'_to_'+str(3-train_round)
            assert int(saved['round_id'])==3-train_round and float(saved['truth'])==r[target]
            assert float(saved['prediction'])==max(0.,float(saved['raw_prediction']))
            prediction_error = max(prediction_error,abs(prediction-float(saved['raw_prediction'])))
            assert np.isclose(prediction,float(saved['raw_prediction']),rtol=1e-8,atol=1e-8)
            checked += 1
    assert checked==len(predictions)==report['prediction_rows'] and len(models)==report['model_fits']==16
    grouped = defaultdict(list)
    for r in predictions:
        grouped[(r['direction'],r['target'],r['strategy'],r['batch_id'])].append(r)
    recomputed = {}
    metric_error = 0.
    for r in read_csv(output/'fold_metrics.csv'):
        key = tuple(r[n] for n in ['direction','target','strategy','batch_id'])
        gs = grouped[key]
        error = np.array([float(p['truth'])-float(p['prediction']) for p in gs])
        m = {'RMSE':float(np.sqrt(np.mean(error**2))),'MAE':float(np.mean(np.abs(error)))}
        assert len(gs)==int(r['points'])
        for name,value in m.items():
            metric_error = max(metric_error,abs(value-float(r[name])))
            assert np.isclose(value,float(r[name]),atol=1e-12,rtol=1e-12)
        recomputed[key] = m
    assert len(recomputed)==len(grouped)==report['fold_scores']==160
    summaries = read_csv(output/'summary.csv')
    assert len(summaries)==len(report['summaries'])==16
    for r,saved in zip(summaries,report['summaries']):
        key = tuple(r[n] for n in ['direction','target','strategy'])
        fs = [m for k,m in recomputed.items() if k[:3]==key]
        ps = [p for p in predictions if tuple(p[n] for n in ['direction','target','strategy'])==key]
        e = np.array([float(p['truth'])-float(p['prediction']) for p in ps])
        expected = {'batch_macro_RMSE':np.mean([f['RMSE'] for f in fs]),'batch_macro_MAE':np.mean([f['MAE'] for f in fs]),
                    'pooled_RMSE':np.sqrt(np.mean(e**2)),'pooled_MAE':np.mean(np.abs(e))}
        assert int(r['batches'])==len(fs)==10 and int(r['pooled_points'])==len(ps)
        assert int(r['negative_raw_predictions'])==sum(float(p['raw_prediction'])<0 for p in ps)
        for name,value in expected.items():
            assert np.isclose(float(r[name]),value,atol=1e-12,rtol=1e-12)
            assert saved[name]==float(r[name])
    pairs = read_csv(output/'paired_batches.csv')
    assert len(pairs)==40
    for r in pairs:
        h = recomputed[(r['direction'],r['target'],'history_ridge',r['batch_id'])]['RMSE']
        i = recomputed[(r['direction'],r['target'],'instantaneous_ridge',r['batch_id'])]['RMSE']
        assert float(r['history_RMSE'])==h and float(r['instantaneous_RMSE'])==i
        assert np.isclose(float(r['delta_RMSE']),h-i,atol=1e-12,rtol=1e-12)
    primary = [r for r in pairs if r['direction']=='1_to_2' and r['target']=='CDW_g_L']
    assert np.isclose(report['primary_comparison']['mean_batch_RMSE_delta'],np.mean([float(r['delta_RMSE']) for r in primary]))
    assert report['primary_comparison']['batches_improved']==sum(float(r['delta_RMSE'])<0 for r in primary)
    result = {'verified':True,'verifier_sha256':sha(__file__),'protocol_sha256':sha(output/'protocol.json'),
              'report_sha256':sha(output/'report.json'),'xlsx_reader':'Independent direct OOXML ZIP/XML parser, no producer import',
              'independent_estimator':'Augmented least squares ridge, independent of producer normal-equation solver',
              'online_rows_independently_parsed':online_points,'offline_rows_independently_parsed':offline_points,
              'all_eligible_features_recomputed':len(independent),'ridge_models_independently_refitted':sum(bool(m['features']) for m in models),
              'prediction_rows_verified':checked,'fold_metrics_verified':len(recomputed),'paired_batches_verified':len(pairs),
              'maximum_feature_absolute_difference':feature_error,'maximum_model_parameter_absolute_difference':model_error,
              'maximum_raw_prediction_absolute_difference':prediction_error,'maximum_fold_metric_absolute_difference':metric_error,
              'actual_sensor_arrival_times_verified':False,'instantaneous_product_concentration_verified':False,'publication_superiority_claim':False}
    (output/'verification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'research/results/ambr_dynamic')
    args = parser.parse_args()
    verify(args.output)
