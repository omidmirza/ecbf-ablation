# Verification of a sample of 50 CSV values against metrics.json

Script `tools/verify_sample.py` (random seed 20261007); exact float equality with the original `final` and `new_anchors` `common/metrics.json` (per-trajectory rows: mean of the 200 CSV values equals the cell mean).

Result: 50 of 50 equal.

| # | file | value checked | CSV | metrics.json | equal |
|---|---|---|---|---|---|
| 1 | summary | final/658ba90110 statistics.h_min.mean | 0.02960064998688655 | 0.02960064998688655 | yes |
| 2 | summary | final/f685599e26 statistics.q_notrans.mean | 0.031848031178969005 | 0.031848031178969005 | yes |
| 3 | summary | new_anchors/c005343555 statistics.m_int.mean | 1.0082245978611255 | 1.0082245978611255 | yes |
| 4 | summary | final/9361a48dcf statistics.q_CBF.mean | 0.05043106945799985 | 0.05043106945799985 | yes |
| 5 | summary | new_anchors/4764ecede6 trajectory_events.qdep.upper95_possible | 0.014867039231272083 | 0.014867039231272083 | yes |
| 6 | summary | final/cfcf266b39 statistics.m_int.mean | 0.6556145893158973 | 0.6556145893158973 | yes |
| 7 | summary | final/b7b67d132b statistics.T_mis.mean | 53.466 | 53.466 | yes |
| 8 | summary | final/5d67b4be57 statistics.m_int.sd | 0.03146029202180308 | 0.03146029202180308 | yes |
| 9 | summary | new_anchors/ff22ce6c7b statistics.q_Y.mean | 0.012287968162876 | 0.012287968162876 | yes |
| 10 | summary | final/03e855ec18 statistics.m_int.sd | 0.08186392843683483 | 0.08186392843683483 | yes |
| 11 | summary | final/dcc08e08b0 statistics.q_use.mean | 0.010461000973838643 | 0.010461000973838643 | yes |
| 12 | summary | final/b2f29f2345 deployment_decision_fraction.possible.mean | 0.0 | 0.0 | yes |
| 13 | summary | new_anchors/637a33b618 statistics.q_unsafe.sd | 0.00027301420123032724 | 0.00027301420123032724 | yes |
| 14 | summary | final/4bcaa3bd66 deployment_decision_fraction.confirmed.mean | 0.0 | 0.0 | yes |
| 15 | summary | final/c9f49143aa statistics.q_nocert.mean | 0.24683047895629448 | 0.24683047895629448 | yes |
| 16 | sweep | T5 corridor phi=0.25 V-C q_CBF sd | 0.012798703778043501 | 0.012798703778043501 | yes |
| 17 | sweep | Default corridor phi=0.5 V-L q_CBF mean | 0.014733119638132428 | 0.014733119638132428 | yes |
| 18 | sweep | Default corridor baseline= c=0 h_min mean | -0.03215387478000547 | -0.03215387478000547 | yes |
| 19 | sweep | Default maze K=50 V-C p_max mean | 0.0469629228211852 | 0.0469629228211852 | yes |
| 20 | sweep | Default corridor phi=0.5 V-U h_min sd | 0.007118060612617034 | 0.007118060612617034 | yes |
| 21 | sweep | Default corridor K=infinity V-C q_dep possible | 0.0 | 0 | yes |
| 22 | sweep | T3 corridor phi=1.0 V-Delta q_use mean | 0.0020658368203654103 | 0.0020658368203654103 | yes |
| 23 | sweep | T2 corridor K=10 V-Delta q_cert possible | 0.0 | 0 | yes |
| 24 | sweep | Default corridor phi=0.75 V-L h_min sd | 0.007232142099486391 | 0.007232142099486391 | yes |
| 25 | sweep | T4 maze delta=0.15 V-L q_dep possible | 0.0 | 0 | yes |
| 26 | sweep | Default corridor phi=0.5 V-C T_mis mean | 35.1365 | 35.1365 | yes |
| 27 | sweep | T3 corridor phi=0.0 V-C T_mis mean | 34.175000000000004 | 34.175000000000004 | yes |
| 28 | sweep | T2 maze phi=0.25 V-Delta m_c sd | 0.07064222455444519 | 0.07064222455444519 | yes |
| 29 | sweep | T2 maze phi=0.5 V-C q_CBF mean | 0.028022962996822765 | 0.028022962996822765 | yes |
| 30 | sweep | T3 maze baseline= ACI q_CBF sd | 0.019575435995645388 | 0.019575435995645388 | yes |
| 31 | heatmap | corridor T_mis a0=0.3 a=0.2 d=0.1 K=50 V-U mean | 34.6225 | 34.6225 | yes |
| 32 | heatmap | corridor q_unsafe a0=0.3 a=0.1 d=0.05 K=10 V-L sd | 0.0 | 0.0 | yes |
| 33 | heatmap | maze m_int a0=0.15 a=0.1 d=0.025 K=200 V-U sd | 0.020671589651926687 | 0.020671589651926687 | yes |
| 34 | heatmap | corridor q_nocert a0=0.3 a=0.1 d=0.15 K=infinity V-C sd | 0.0593340301851933 | 0.0593340301851933 | yes |
| 35 | heatmap | corridor runtime_ms a0=0.3 a=0.05 d=0.2 K=50 V-L mean | 1.4537257941567798 | 1.4537257941567798 | yes |
| 36 | heatmap | corridor q_notrans a0=0.1 a=0.05 d=0.05 K=infinity V-Delta sd | 0.2404364414682821 | 0.2404364414682821 | yes |
| 37 | heatmap | maze m_Y a0=0.25 a=0.1 d=0.15 K=50 V-L mean | 0.16991471053412002 | 0.16991471053412002 | yes |
| 38 | heatmap | corridor q_CBF a0=0.3 a=0.1 d=0.1 K=100 V-Delta sd | 0.009978473478506337 | 0.009978473478506337 | yes |
| 39 | heatmap | corridor q_CBF a0=0.25 a=0.1 d=0.025 K=50 V-C sd | 0.010708236412093954 | 0.010708236412093954 | yes |
| 40 | heatmap | maze m_int a0=0.25 a=0.15 d=0.05 K=100 V-U sd | 0.025321469387354092 | 0.025321469387354092 | yes |
| 41 | heatmap | maze m_c a0=0.2 a=0.15 d=0.05 K=100 V-Delta sd | 0.0 | 0.0 | yes |
| 42 | heatmap | maze m_c a0=0.2 a=0.1 d=0.05 K=100 V-C mean | 0.8609222515559417 | 0.8609222515559417 | yes |
| 43 | heatmap | maze p_max a0=0.25 a=0.1 d=0.05 K=infinity V-L mean | 0.04636831060628611 | 0.04636831060628611 | yes |
| 44 | heatmap | corridor runtime_ms a0=0.2 a=0.15 d=0.05 K=10 V-U mean | 1.7323409620557595 | 1.7323409620557595 | yes |
| 45 | heatmap | maze q_nocert a0=0.25 a=0.15 d=0.15 K=infinity V-L mean | 0.11224042605906794 | 0.11224042605906794 | yes |
| 46 | per-trajectory | T4 corridor ACI T_mis: mean of 200 CSV values | 33.128 | 33.128 | yes |
| 47 | per-trajectory | T1 corridor ACI q_CBF: mean of 200 CSV values | 0.17675774490778098 | 0.17675774490778098 | yes |
| 48 | per-trajectory | T5 maze Zero margin q_CBF: mean of 200 CSV values | 0.24934838971844625 | 0.24934838971844625 | yes |
| 49 | per-trajectory | Default corridor V-U h_min: mean of 200 CSV values | 0.048896080303046534 | 0.048896080303046534 | yes |
| 50 | per-trajectory | T1 maze ACI q_CBF: mean of 200 CSV values | 0.1583983083953066 | 0.1583983083953066 | yes |
