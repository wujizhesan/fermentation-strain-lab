import hashlib
import unittest

import numpy as np

from ambr_dynamic import CHANNELS, RECIPE, feature_names, fit, metadata, partitions, predict, snapshot, validate_source


class AmbrDynamicTests(unittest.TestCase):
    def signal(self):
        return np.array([[i+2,float(i)]+[float(i+j) for j in range(len(CHANNELS))] for i in range(8)])

    def row(self, time):
        return {**metadata('T2 YE 0.5 (1)'),'time_h':time,'row_id':'x','CDW_g_L':2.}

    def test_future_signal_mutation_and_extension_do_not_change_snapshot(self):
        a = self.signal()
        b = a.copy()
        b[b[:,1]>3.25,2:] = 1e9
        b = np.vstack([b,np.array([99,99]+[1e10]*len(CHANNELS))])
        first,_ = snapshot(a,self.row(3.25))
        second,_ = snapshot(b,self.row(3.25))
        self.assertEqual(first,second)
        self.assertEqual(first['prefix_points'],4)

    def test_exact_boundary_and_no_future_interpolation(self):
        a = self.signal()
        first,_ = snapshot(a,self.row(3.))
        self.assertEqual(first['last_CER'],3.)
        second,_ = snapshot(a,self.row(2.999999999),max_age_h=1.)
        self.assertEqual(second['last_CER'],2.)
        self.assertEqual(second['mean_CER_1h'],2.)

    def test_empty_prefix_and_stale_signal_are_rejected(self):
        self.assertEqual(snapshot(self.signal(),self.row(-.0-.0))[-1],'')
        first,reason = snapshot(self.signal()[1:],self.row(.5))
        self.assertIsNone(first)
        self.assertEqual(reason,'no_past_online_observation')
        self.assertIsNone(snapshot(self.signal(),self.row(8.)) [0])

    def test_snapshot_ignores_assays_and_offline_final_volume(self):
        r = self.row(3.)
        expected = snapshot(self.signal(),r)
        r.update({'CDW_g_L':10000.,'mucic_acid_corrected_g_L':999999.,'offline_final_volume':1e20})
        self.assertEqual(snapshot(self.signal(),r),expected)
        self.assertFalse(any('acid' in n or 'base' in n or n=='last_volume' for n in feature_names('history_ridge')))

    def test_round_holdout_preserves_cultures_and_missing_target(self):
        rows = [{**self.row(2.),'row_id':str(i),'round_id':1 if i<2 else 2,'batch_id':'a' if i<2 else 'b',
                 'feature_eligible':True,'CDW_g_L':None if i==1 else 2.} for i in range(4)]
        train,test = partitions(rows,1,'CDW_g_L')
        self.assertEqual([r['row_id'] for r in train],['0'])
        self.assertEqual([r['row_id'] for r in test],['2','3'])
        with self.assertRaises(ValueError):
            partitions([{**r,'batch_id':'shared'} for r in rows],1,'CDW_g_L')

    def test_training_normalization_and_query_label_independence(self):
        rows = [{**self.row(float(i)),**{n:float(i) for n in RECIPE},'row_id':str(i),'CDW_g_L':2.+3*i} for i in range(5)]
        model = fit(rows,'CDW_g_L','time_recipe_ridge',10.)
        self.assertEqual(model['center'],[2.]*6)
        query = {**rows[2],'CDW_g_L':-1e10,'round_id':2}
        np.testing.assert_array_equal(predict(model,[rows[2]]),predict(model,[query]))
        self.assertEqual(model,fit(rows,'CDW_g_L','time_recipe_ridge',10.))

    def test_replicate_identity_is_explicit_and_unknown_not_inferred(self):
        r = metadata('D-161646 pH 3.5 (2)')
        self.assertEqual((r['round_id'],r['strain'],r['pH_setpoint']),(2,'D-161646',3.5))
        r = metadata('T2 glucose 10 (1)')
        self.assertEqual((r['YE_nominal_g_L'],r['glucose_feed_nominal_g_L']),(1.,10.))
        with self.assertRaises(ValueError):
            metadata('unknown culture (1)')

    def test_source_tampering_is_rejected(self):
        data = b'original'
        spec = {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'published_md5':hashlib.md5(data).hexdigest()}
        validate_source(data,spec)
        with self.assertRaises(ValueError):
            validate_source(b'changed',spec)

    def test_integrals_and_window_means_use_only_recorded_prefix(self):
        features,_ = snapshot(self.signal(),self.row(3.25))
        self.assertEqual(features['integral_CER'],4.5)
        self.assertEqual(features['mean_CER_6h'],1.5)
        self.assertEqual(features['mean_CER_1h'],3.)
        with self.assertRaises(ValueError):
            snapshot(self.signal()[::-1],self.row(3.))


if __name__=='__main__':
    unittest.main()
