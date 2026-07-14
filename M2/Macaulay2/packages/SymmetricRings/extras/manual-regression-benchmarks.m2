-- Manual cold regression benchmarks for SymmetricRings.
--
-- This file is intentionally separate from the systematic suite. Each benchmark is a
-- self-contained block: select the whole block and run it.  The SymmetricRings
-- blocks use frac(QQ[t]) by default, because that is the coefficient field
-- most likely to expose regressions in the current code.
--
-- Timings below were taken with the local Codex M2 build on 2026-07-05.
-- Comments also record earlier cold QQ timings for SymmetricRings and cold QQ
-- timings for SchurRings when there is an analogous Schur-in/Schur-out test.
-- A timing such as ">60s" means the cold block reached that timeout cap in
-- the rerun, not that the exact runtime is known.
--
-- When rerunning the benchmarks, DO NOT overwrite the previous timings.
-- Add a new comment, e.g.:
-- "-- 2026-07-15 04:32pm; codex frac(QQ[t]): 0.841233s; codex QQ: 0.761376s; SchurRings QQ: 9.51899s"
-- If there are multiple "time" tests in a block, add a separate comment after each "time"

------------------------
-- Plethym benchmarks --
------------------------

-- 01. Schur plethysm, two-row outer and inner.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{5,1}@S_{3,1};             -- codex frac(QQ[t]): 0.841233s; codex QQ: 0.761376s; SchurRings QQ: 9.51899s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.70741s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 1.10895s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 1.37243s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 1.17262s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.824578s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.760727s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.755312s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.781392s; 4-run median: 0.78048s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.759495s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.759635s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.75614s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.776925s; runs: 0.776113s, 0.776925s, 0.796664s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{5,1}@s_{3,1};             -- codex SchurRings QQ: 9.51899s

-- 02. Same plethysm split into p-plethysm and conversion.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time F5131 = plethysm(S_{5,1}, S_{3,1}); -- codex frac(QQ[t]): 0.006968s; codex QQ: 0.001481s; SchurRings QQ plethysm: 9.32736s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.007569s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.007586s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.008676s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.007837s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.007604s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.007391s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.007392s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.00757s; 4-run median: 0.0074785s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.006878s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.0071s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.006793s
-- 2026-07-12; 3-run median: 0.007201s; runs: 0.007347s, 0.007201s, 0.007164s
time toS F5131;                   -- codex frac(QQ[t]): 1.45196s; codex QQ: 0.765673s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.723232s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 2.73161s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 3.43703s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 2.80095s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 2.82065s
-- 2026-07-10 12:43pm EDT; codex frac(QQ[t]): 1.38453s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 1.32541s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 1.33186s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 1.40884s; 4-run median: 1.322515s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 1.33755s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 1.36907s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.754161s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.745158s; runs: 0.750009s, 0.745158s, 0.740663s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time plethysm(s_{5,1}, s_{3,1});  -- codex SchurRings QQ: 9.32736s

-- 03. One-row plethysm with expensive Schur conversion.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_8@S_3;                     -- codex frac(QQ[t]): 2.43985s; codex QQ: 2.42998s; SchurRings QQ: 8.94662s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.810654s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 2.4508s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 3.13797s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 2.45078s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 2.46251s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.889134s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.886118s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.870545s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.898218s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.88368s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.870136s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 3-run median: 0.823794s; runs: 0.792454s, 0.823794s, 0.832865s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 4.21643s; runs: 4.21643s, 4.20923s, 4.23902s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_8@s_3;                     -- codex SchurRings QQ: 8.94662s

-- 04. Same one-row plethysm split into p-plethysm and conversion.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time F83 = plethysm(S_8, S_3);    -- codex frac(QQ[t]): 0.015178s; codex QQ: 0.00321s; SchurRings QQ plethysm: 9.06517s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.015899s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.01521s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.018059s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.016171s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.015939s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.015645s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.01527s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.017048s; 4-run median: 0.015422s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.015161s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.016591s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.015072s
-- 2026-07-12; 3-run median: 0.015232s; runs: 0.015232s, 0.015433s, 0.015029s
time toS F83;                     -- codex frac(QQ[t]): 2.87714s; codex QQ: 0.981811s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.806766s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 9.77278s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 12.6455s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 10.3423s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 10.7248s
-- 2026-07-10 12:43pm EDT; codex frac(QQ[t]): 2.48308s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 2.49381s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 2.4381s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 2.45434s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 2.59383s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 2.63768s; 4-run median: 2.617725s
-- 2026-07-10 4:16pm EDT; codex frac(QQ[t]): 0.876586s; 4-run median: 0.8871075s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.90094s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 3.97126s; runs: 4.26755s, 3.97126s, 3.94282s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time plethysm(s_8, s_3);          -- codex SchurRings QQ: 9.06517s

-- 05. Four-part outer partition against S_2.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{6,3,2,1}@S_2;             -- codex frac(QQ[t]): 1.10539s; codex QQ: 1.06687s; SchurRings QQ: 10.537s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 1.08674s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 1.14527s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 1.46592s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 1.15597s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 1.15337s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 1.09818s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 1.08479s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 1.09346s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 1.08895s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 1.10869s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 1.09283s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 1.09263s; runs: 1.08904s, 1.09263s, 1.0946s
-- 2026-07-12 6:30PM EDT; shape-aware post-plethysm 3-run median: 0.752006s; runs: 0.74745s, 0.752006s, 0.754022s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{6,3,2,1}@s_2;             -- codex SchurRings QQ: 10.537s

-- 06. Larger four-part outer partition against S_2.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{5,4,3,1}@S_2;             -- codex frac(QQ[t]): 3.97659s; codex QQ: 3.8739s; SchurRings QQ: 72.6661s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 3.94761s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 4.08913s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 5.20889s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 4.05939s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 4.06262s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 3.98741s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 3.96723s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 3.9724s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 3.96322s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 3.96801s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 3.98746s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 3.97757s; runs: 3.97757s, 3.96596s, 4.02501s
-- 2026-07-12 6:30PM EDT; shape-aware post-plethysm 3-run median: 2.11634s; runs: 2.0962s, 2.11634s, 2.14892s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{5,4,3,1}@s_2;             -- codex SchurRings QQ: 72.6661s

-- 07. Midrange S_2 plethysm.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{7,4,2}@S_2;               -- codex frac(QQ[t]): 0.452516s; codex QQ: 0.342234s; SchurRings QQ: 68.1274s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.353646s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.491661s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.597604s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.4982s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.492872s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.359331s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.351816s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.353176s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.354669s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.358s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.356501s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.357914s; runs: 0.362777s, 0.357914s, 0.357498s
-- 2026-07-12 6:30PM EDT; retained Adams/Jacobi-Trudi 3-run median: 0.360677s; runs: 0.358662s, 0.360677s, 0.371452s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{7,4,2}@s_2;               -- codex SchurRings QQ: 68.1274s

-- 08. Larger S_2 plethysm.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{8,4,2}@S_2;               -- codex frac(QQ[t]): 0.911587s; codex QQ: 0.700197s; SchurRings QQ: >60s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.744928s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.923382s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 1.23764s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.936613s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.936267s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.760856s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.73614s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.74835s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.715677s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.734768s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.739211s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.772101s; runs: 0.769441s, 0.772101s, 0.773263s
-- 2026-07-12 6:30PM EDT; retained Adams/Jacobi-Trudi 3-run median: 0.759063s; runs: 0.754231s, 0.759063s, 0.76743s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{8,4,2}@s_2;               -- codex SchurRings QQ: >60s

-- 09. Another S_2 plethysm near one second over frac(QQ[t]).
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{7,4,3}@S_2;               -- codex frac(QQ[t]): 0.932549s; codex QQ: 0.7003s; SchurRings QQ: >60s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.74959s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.929912s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 1.11841s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.947592s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.939982s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.778702s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.74582s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.746429s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.723545s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.727437s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.725523s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.765776s; runs: 0.765776s, 0.763207s, 0.765859s
-- 2026-07-12 6:30PM EDT; retained Adams/Jacobi-Trudi 3-run median: 0.755319s; runs: 0.755028s, 0.755319s, 0.760999s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{7,4,3}@s_2;               -- codex SchurRings QQ: >60s

-- 10. S_{2,1} plethysm, moderate.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{5,3}@S_{2,1};             -- codex frac(QQ[t]): 0.836392s; codex QQ: 0.744417s; SchurRings QQ: 12.4018s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.735081s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.868491s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 1.10282s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.894284s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.889246s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.766393s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.730682s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.744864s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.749526s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.738693s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.730341s
-- 2026-07-11 2:00pm EDT; codex frac(QQ[t]) 3-run median: 0.710339s; runs: 0.705638s, 0.710339s, 0.73303s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.494502s; runs: 0.496404s, 0.494502s, 0.471265s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{5,3}@s_{2,1};             -- codex SchurRings QQ: 12.4018s

-- 11. S_{2,1} plethysm, smaller but still useful.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{4,3}@S_{2,1};             -- codex frac(QQ[t]): 0.233528s; codex QQ: 0.214275s; SchurRings QQ: 6.26026s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.24138s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.254874s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.321636s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.254952s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.26207s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.21513s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.21289s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.216756s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.216524s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.214346s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.220832s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.121509s; runs: 0.121156s, 0.124889s, 0.121509s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{4,3}@s_{2,1};             -- codex SchurRings QQ: 6.26026s

-- 12. S_{3,1} plethysm, smaller reference case.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{4,2}@S_{3,1};             -- codex frac(QQ[t]): 0.686352s; codex QQ: 0.595111s; SchurRings QQ: 7.08531s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.578886s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.734983s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.942472s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.743213s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.749864s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.589368s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.592511s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.593199s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.596664s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.605893s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.588963s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 0.456437s; runs: 0.444452s, 0.456437s, 0.471189s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{4,2}@s_{3,1};             -- codex SchurRings QQ: 7.08531s

-- 13. S_{2,1} plethysm, heavier.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{5,4}@S_{2,1};             -- codex frac(QQ[t]): 2.53696s; codex QQ: 2.29245s; SchurRings QQ: 45.8267s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 2.28557s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 2.76183s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 3.28681s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 2.77302s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 2.50091s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 2.32009s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 2.31306s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 2.31188s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 2.31199s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 2.28687s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 2.28301s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 1.70851s; runs: 1.70851s, 1.6258s, 1.74626s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{5,4}@s_{2,1};             -- codex SchurRings QQ: 45.8267s

-- 14. Three-part outer partition against S_{2,1}.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{5,3,1}@S_{2,1};           -- codex frac(QQ[t]): 2.52421s; codex QQ: 2.28816s; SchurRings QQ: 30.0077s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 2.00866s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 2.1918s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 2.67296s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 2.24947s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 2.2472s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 2.05917s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 2.01777s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 2.0045s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 2.03933s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 2.06803s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 2.08132s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 1.22687s; runs: 1.22687s, 1.21875s, 1.23168s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{5,3,1}@s_{2,1};           -- codex SchurRings QQ: 30.0077s

-- 15. Another three-part outer partition against S_{2,1}.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{4,3,2}@S_{2,1};           -- codex frac(QQ[t]): 2.64538s; codex QQ: 2.39574s; SchurRings QQ: 53.2158s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 2.37581s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 2.87853s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 3.46437s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 3.01785s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 2.6182s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 2.4365s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 2.41405s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 2.41226s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 2.42425s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 2.39511s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 2.39687s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 1.86705s; runs: 1.86705s, 1.87605s, 1.84993s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{4,3,2}@s_{2,1};           -- codex SchurRings QQ: 53.2158s

-- 16. S_{3,1} plethysm, heavier.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{5,2}@S_{3,1};             -- codex frac(QQ[t]): 3.70285s; codex QQ: 3.42277s; SchurRings QQ: 30.7652s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 3.2399s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 6.87293s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 8.1419s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 7.12635s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 3.59857s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 3.31066s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 3.32447s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 3.31868s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 3.31196s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 3.3194s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 3.31126s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 5.82679s; runs: 5.84026s, 5.82679s, 5.57444s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{5,2}@s_{3,1};             -- codex SchurRings QQ: 30.7652s

-- 17. S_{3,1} plethysm, heavier.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time S_{4,3}@S_{3,1};             -- codex frac(QQ[t]): 4.09575s; codex QQ: 3.78781s; SchurRings QQ: 39.3563s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 3.63223s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 7.70281s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 9.2059s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 7.9734s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 3.95101s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 3.7374s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 3.67974s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 3.68167s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 3.70733s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 3.6876s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 3.67192s
-- 2026-07-12; automatic p-to-S dispatcher 3-run median: 6.22708s; runs: 6.06079s, 6.22708s, 6.25281s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{4,3}@s_{3,1};             -- codex SchurRings QQ: 39.3563s

-- 18. Power-sum to Schur conversion.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toS p_{9,7,5,3};             -- codex frac(QQ[t]): 0.52265s; codex QQ: 0.500656s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.034835s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.079701s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.098224s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.034765s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.034575s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.035137s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.033894s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.033998s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.034479s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.03427s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.033248s
-- 2026-07-11 2:00pm EDT; codex frac(QQ[t]) 3-run median: 0.033542s; runs: 0.032698s, 0.034564s, 0.033542s
-- 2026-07-12; automatic abacus route 3-run median: 0.004616s; runs: 0.004616s, 0.004614s, 0.004763s

-- 19. Power-sum plethysm followed by Schur conversion.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toS(p_{7,5,3} @ p_2);        -- codex frac(QQ[t]): 3.63397s; codex QQ: 3.61262s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.126104s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.173517s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.219411s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.124512s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.124869s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.12991s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.125915s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.127086s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.127757s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.12777s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.125829s
-- 2026-07-12; automatic post-plethysm abacus route 3-run median: 0.006941s; runs: 0.006772s, 0.006941s, 0.007114s

-- 20. Another power-sum plethysm followed by Schur conversion.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toS(p_{5,4,3,2} @ p_2);      -- codex frac(QQ[t]): 2.08908s; codex QQ: 1.94596s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.114543s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.219157s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.261403s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.112172s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.112323s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.116508s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.113887s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.114641s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.114775s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.113051s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.113328s
-- 2026-07-12; automatic post-plethysm abacus route 3-run median: 0.01152s; runs: 0.01152s, 0.011379s, 0.011699s

--------------------
-- Schur products --
--------------------

-- 21. Schur product followed by toS; the SchurRings analogue is just product.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toS(S_{12,8}*S_{10,6});      -- codex frac(QQ[t]): 0.001705s; codex QQ: 0.001184s; SchurRings QQ product: 0.000311s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.001598s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.001928s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.009548s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.001449s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.001779s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.00146s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.001425s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.001403s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.001363s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.001431s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.001423s
-- 2026-07-11 2:00pm EDT; codex frac(QQ[t]) 3-run median: 0.001339s; runs: 0.00133s, 0.001339s, 0.001388s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{12,8}*s_{10,6};           -- codex SchurRings QQ: 0.000311s

-- 22. Larger Schur product followed by toS.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toS(S_{10,8,4}*S_{9,7,3});   -- codex frac(QQ[t]): 0.014058s; codex QQ: 0.009821s; SchurRings QQ product: 0.009201s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.014015s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.014825s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.167704s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.013645s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.012812s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.012712s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.012523s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.01217s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.011706s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.01154s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.011771s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{10,8,4}*s_{9,7,3};        -- codex SchurRings QQ: 0.009201s

-- 23. Four-part Schur product followed by toS.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toS(S_{8,6,4,2}*S_{7,5,3,1}); -- codex frac(QQ[t]): 0.02398s; codex QQ: 0.018298s; SchurRings QQ product: 0.028526s
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]) 3-run median: 0.021327s; runs: 0.02647s, 0.020424s, 0.021327s
-- 2026-07-09 10:56pm EDT; codex frac(QQ[t]): 0.024165s
-- 2026-07-10 11:28pm EDT; codex frac(QQ[t]): 0.318697s
-- 2026-07-10 11:49pm EDT; codex frac(QQ[t]): 0.022052s
-- 2026-07-10 12:11pm EDT; codex frac(QQ[t]): 0.021733s
-- 2026-07-10 1:42pm EDT; codex frac(QQ[t]): 0.021369s
-- 2026-07-10 2:21pm EDT; codex frac(QQ[t]): 0.020891s
-- 2026-07-10 2:49pm EDT; codex frac(QQ[t]): 0.020799s
-- 2026-07-10 3:37pm EDT; codex frac(QQ[t]): 0.01999s
-- 2026-07-10 4:07pm EDT; codex frac(QQ[t]): 0.020169s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.020884s

restart
needsPackage "SchurRings"
SR = schurRing(QQ, s, 40)
time s_{8,6,4,2}*s_{7,5,3,1};     -- codex SchurRings QQ: 0.028526s

-- 24. Retained Hall-Littlewood product through multiplicative generators.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time multiplyToBasis(Q_{4,2}, Q_{3,1}, Q);
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.011132s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.01408s
-- 2026-07-11 2:00pm EDT; codex frac(QQ[t]) 3-run median: 0.008984s; runs: 0.008617s, 0.008984s, 0.008985s

-- 25. The same Hall-Littlewood product after operand grouping is lost.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(Q_{4,2}*Q_{3,1}, Q);
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.068011s
-- 2026-07-10 4:44pm EDT; codex frac(QQ[t]): 0.11608s

---------------------------
-- Plethysm coefficients --
---------------------------

-- 26. Research benchmark: Hall-Littlewood plethysm paired with a normalized
-- Hall-Littlewood function of the same weight.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time hallInnerProduct(plethysm(Q_{1,1}, Q_{8,2}), P_{14,4,2})
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 1.88454s
-- 2026-07-10 7:08pm EDT; codex frac(QQ[t]) 3-run median: 1.97601s; runs: 1.98776s, 1.89457s, 1.97601s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 1.89875s; runs: 1.89875s, 1.90417s, 1.84039s
-- 2026-07-11 8:59am EDT; codex frac(QQ[t]) 3-run median: 0.776079s; runs: 0.776079s, 0.775621s, 0.785946s
-- 2026-07-12 6:30PM EDT; targeted dual coefficient 3-run median: 0.729709s; runs: 0.708903s, 0.729709s, 0.758376s

-- 27. Lower-weight version of the research benchmark.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time hallInnerProduct(plethysm(Q_{1,1}, Q_{6,2}), P_{10,4,2})
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.511706s
-- 2026-07-10 7:08pm EDT; codex frac(QQ[t]) 3-run median: 0.544034s; runs: 0.550225s, 0.538848s, 0.544034s
-- 2026-07-11 8:59am EDT; codex frac(QQ[t]) 3-run median: 0.188873s; runs: 0.190829s, 0.188873s, 0.187281s

-- 28. One-row outer Hall-Littlewood plethysm and Hall inner product.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time hallInnerProduct(plethysm(Q_2, Q_{6,2}), P_{12,4})
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 0.411301s
-- 2026-07-10 7:08pm EDT; codex frac(QQ[t]) 3-run median: 0.41303s; runs: 0.411968s, 0.41303s, 0.419661s
-- 2026-07-11 8:59am EDT; codex frac(QQ[t]) 3-run median: 0.143973s; runs: 0.135147s, 0.159811s, 0.143973s

-- 29. Three-box outer Hall-Littlewood plethysm and Hall inner product.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time hallInnerProduct(plethysm(Q_{2,1}, Q_{4,2}), P_{10,5,3})
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]): 1.03807s
-- 2026-07-10 7:08pm EDT; codex frac(QQ[t]) 3-run median: 1.13079s; runs: 1.12368s, 1.1498s, 1.13079s
-- 2026-07-11 8:59am EDT; codex frac(QQ[t]) 3-run median: 0.414548s; runs: 0.414548s, 0.416405s, 0.409409s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 3-run median: 0.381229s; runs: 0.377035s, 0.381229s, 0.382835s

-- 30. Omega-dual Hall-Littlewood plethysm and Hall inner product.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time hallInnerProduct(plethysm(B_{1,1}, B_{8,2}), Pomega_{14,4,2})
-- 2026-07-10 7:54pm EDT; codex frac(QQ[t]) 3-run median: 1.89154s; runs: 1.89154s, 1.87699s, 1.90622s
-- 2026-07-10 7:08pm EDT; codex frac(QQ[t]) 3-run median: 0.367093s; runs: 0.367093s, 0.38069s, 0.365832s
-- 2026-07-11 8:59am EDT; codex frac(QQ[t]) 3-run median: 0.783413s; runs: 0.783413s, 0.775347s, 0.805431s
-- 2026-07-12 6:30PM EDT; targeted dual coefficient 3-run median: 0.738145s; runs: 0.708887s, 0.738145s, 0.757765s

----------------------------------------------------------------------------------
-- Built-in source-basis-to-power-sum benchmarks. The power-sum identity route  --
-- is omitted; each nontrivial conversion is sized for a roughly 0.2s-1.0s run. --
----------------------------------------------------------------------------------

-- 31. Complete homogeneous to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(h_30, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.486405s; runs: 0.493359s, 0.486405s, 0.478583s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.097663s; runs: 0.111517s, 0.097467s, 0.097663s
-- 2026-07-11 2:00pm EDT; codex frac(QQ[t]) 3-run median: 0.070735s; runs: 0.070883s, 0.070735s, 0.068267s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 6-run median: 0.0941715s; runs: 0.090954s, 0.095183s, 0.092969s, 0.096455s, 0.101884s, 0.09316s
-- 2026-07-12 6:30PM EDT; automatic QQ shadow 3-run median: 0.074349s; runs: 0.0741s, 0.074349s, 0.075743s

-- 32. Elementary to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(e_30, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.467567s; runs: 0.460401s, 0.477532s, 0.467567s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.098054s; runs: 0.098054s, 0.099104s, 0.095947s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 3-run median: 0.10011s; runs: 0.110485s, 0.099881s, 0.10011s
-- 2026-07-12 6:30PM EDT; automatic QQ shadow 3-run median: 0.07436s; runs: 0.071265s, 0.07436s, 0.086908s

-- 33. Monomial expansion to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(m_{5,3,2} + m_{6,3,1}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.270322s; runs: 0.267037s, 0.272572s, 0.270322s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.000411s; runs: 0.000412s, 0.000405s, 0.000411s

-- 34. Forgotten expansion to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(ff_{5,3,2} + ff_{6,3,1}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.27273s; runs: 0.269645s, 0.27273s, 0.276457s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.000435s; runs: 0.000435s, 0.000438s, 0.000399s

-- 35. Schur to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(S_{13,9,6}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.319999s; runs: 0.324172s, 0.319999s, 0.317388s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.198488s; runs: 0.196919s, 0.199921s, 0.198488s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 3-run median: 0.194491s; runs: 0.193336s, 0.19544s, 0.194491s

-- 36. Default-normalized Schur Omega construction to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(Somega_{11,8,6}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.688168s; runs: 0.683136s, 0.691278s, 0.688168s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.65738s; runs: 0.664829s, 0.65738s, 0.65305s

-- 37. Hall-Littlewood q generator to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(q_21, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.260543s; runs: 0.260543s, 0.260968s, 0.257525s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.213096s; runs: 0.210341s, 0.213096s, 0.21389s

-- 38. Hall-Littlewood b generator to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(b_21, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.263863s; runs: 0.265672s, 0.260861s, 0.263863s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.211013s; runs: 0.212646s, 0.204965s, 0.211013s

-- 39. Capital Hall-Littlewood Q to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(Q_{14,4,2}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.503719s; runs: 0.503719s, 0.502475s, 0.504028s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.406321s; runs: 0.402219s, 0.406321s, 0.406917s

-- 40. Capital Hall-Littlewood B to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(B_{14,4,2}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.504564s; runs: 0.513946s, 0.504043s, 0.504564s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.405433s; runs: 0.404738s, 0.413585s, 0.405433s

-- 41. Normalized Hall-Littlewood P to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(P_{14,4,2}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.614821s; runs: 0.623244s, 0.612292s, 0.614821s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.536961s; runs: 0.526327s, 0.557115s, 0.536961s

-- 42. Hall-Littlewood P Omega to power sums.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(Pomega_{14,4,2}, p)
-- 2026-07-10 7:15pm EDT; codex frac(QQ[t]) 3-run median: 0.624045s; runs: 0.625423s, 0.624045s, 0.61519s
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.538536s; runs: 0.541402s, 0.524746s, 0.538536s
-- 2026-07-11 2:00pm EDT; codex frac(QQ[t]) 3-run median: 0.542237s; runs: 0.534633s, 0.546879s, 0.542237s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 3-run median: 0.530273s; runs: 0.541525s, 0.530273s, 0.522459s

-- 43. Unnormalized Schur Omega to power sums through the actual Somega kernel.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]), "NormalizeSomega" => false)
time toBasis(Somega_{14,10,6}, p)
-- 2026-07-10 7:26pm EDT; codex frac(QQ[t]) 3-run median: 0.317636s; runs: 0.314168s, 0.317636s, 0.318478s

-- Built-in power-sum-to-target-basis benchmarks. The power-sum identity route
-- is omitted. The dense p -> m/ff cases sit just below the intended timing
-- range because increasing their length by one raises the runtime above 1s.

-- 44. Power sums to complete homogeneous.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_30, h)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.379017s; runs: 0.395746s, 0.379017s, 0.376827s
-- 2026-07-10 7:45pm EDT; pre-GMP overflow fix, invalid result timing; 3-run median: 0.075551s; runs: 0.073186s, 0.075551s, 0.07647s
-- 2026-07-10 7:48pm EDT; codex frac(QQ[t]) 3-run median: 0.079341s; runs: 0.094787s, 0.07612s, 0.079341s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 3-run median: 0.075536s; runs: 0.089334s, 0.074708s, 0.075536s

-- 45. Power sums to elementary.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_30, e)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.492337s; runs: 0.492337s, 0.496361s, 0.488107s
-- 2026-07-10 7:45pm EDT; pre-GMP overflow fix, invalid result timing; 3-run median: 0.075546s; runs: 0.075546s, 0.075933s, 0.074591s
-- 2026-07-10 7:48pm EDT; codex frac(QQ[t]) 3-run median: 0.07665s; runs: 0.080989s, 0.07665s, 0.076064s

-- 46. Dense power sum to monomial.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_{1,1,1,1,1,1,1,1,1,1}, m)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.117548s; runs: 0.116233s, 0.117548s, 0.119597s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.000897s; runs: 0.000933s, 0.000887s, 0.000897s

-- 47. Dense power sum to forgotten.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_{1,1,1,1,1,1,1,1,1,1}, ff)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.117639s; runs: 0.117281s, 0.117639s, 0.12122s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.00089s; runs: 0.000808s, 0.000956s, 0.00089s

-- 48. Power sums to Schur.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_{13,8,6,3}, S)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.235515s; runs: 0.234832s, 0.235515s, 0.236179s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.236601s; runs: 0.236601s, 0.237423s, 0.236213s
-- 2026-07-11 2:00pm EDT; codex frac(QQ[t]) 3-run median: 0.230702s; runs: 0.229068s, 0.230702s, 0.231927s
-- 2026-07-12 12:28pm EDT; codex frac(QQ[t]) 3-run median: 0.234s; runs: 0.24197s, 0.232648s, 0.234s
-- 2026-07-12; automatic abacus route 3-run median: 0.010708s
-- 2026-07-12; automatic abacus route 3-run median: 0.010516s; runs: 0.010516s, 0.010527s, 0.010384s

-- 48a. The same power-sum-to-Schur input explicitly forced to use beta-set/
-- abacus rim-hook generation. Run this block in a fresh process with:
-- M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=abacus-rim-hooks M2
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_{13,8,6,3}, S)
-- 2026-07-12; codex frac(QQ[t]) first comparison: generate/filter 0.241353s; abacus 0.010526s
-- 2026-07-12; fresh-process 3-run medians: generate/filter 0.230751s; abacus 0.010708s
-- 2026-07-12; forced abacus route 3-run median: 0.010156s; runs: 0.011096s, 0.010123s, 0.010156s

-- Correctness sweep for the forced abacus route. This compares every
-- power-sum basis element through weight 12 with the independent p -> h -> S
-- route in the same process.
scan(0..12, n -> scan(partitions n, mu -> (
            F := p_(toList mu);
            assert(toS F == toBasis(toBasis(F, h), S))
            )))

-- 49. Power sums to the unnormalized Schur Omega basis.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]), "NormalizeSomega" => false)
time toBasis(p_{9,7,5,3}, Somega)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.33743s; runs: 0.335309s, 0.33743s, 0.337661s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.032346s; runs: 0.032463s, 0.032346s, 0.032301s
-- 2026-07-12; automatic omega-abacus route 3-run median: 0.00278s; runs: 0.00278s, 0.002712s, 0.002815s

-- 50. Power sums to Hall-Littlewood q generators.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_18, q)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.61404s; runs: 0.612093s, 0.61404s, 0.618211s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.007403s; runs: 0.007398s, 0.007403s, 0.007579s

-- 51. Power sums to Hall-Littlewood b generators.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_18, b)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.613664s; runs: 0.610264s, 0.613664s, 0.620897s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.007369s; runs: 0.007424s, 0.007278s, 0.007369s

-- 52. Power sums to capital Hall-Littlewood Q.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_8, Q)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.376881s; runs: 0.379144s, 0.376881s, 0.367133s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.295576s; runs: 0.295647s, 0.295576s, 0.289414s

-- 53. Power sums to capital Hall-Littlewood B.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_8, B)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.36959s; runs: 0.36959s, 0.366984s, 0.377864s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.299394s; runs: 0.301326s, 0.298614s, 0.299394s

-- 54. Power sums to normalized Hall-Littlewood P.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_8, P)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.378564s; runs: 0.373857s, 0.378564s, 0.385915s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.298312s; runs: 0.289761s, 0.299562s, 0.298312s

-- 55. Power sums to Hall-Littlewood P Omega.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time toBasis(p_8, Pomega)
-- 2026-07-10 7:32pm EDT; codex frac(QQ[t]) 3-run median: 0.373282s; runs: 0.37098s, 0.373282s, 0.374472s
-- 2026-07-10 7:45pm EDT; codex frac(QQ[t]) 3-run median: 0.2966s; runs: 0.299602s, 0.2966s, 0.28474s

----------------------------------------------------------------------------
-- Comparative inner-product pipeline benchmarks. Run a block in a fresh M2
-- process for a cold cache. The second call records the corresponding warm
-- cache, and reversed calls check symmetric route selection.
----------------------------------------------------------------------------

-- 56. Schur-complete Kostka route, cold/warm and both orientations.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing QQ
time hallInnerProduct(S_{8,5,3}, h_{5,4,3,2,1,1})
time hallInnerProduct(S_{8,5,3}, h_{5,4,3,2,1,1})
time hallInnerProduct(h_{5,4,3,2,1,1}, S_{8,5,3})
-- 2026-07-11; codex QQ cold/warm/reversed: 0.001214s, 0.000182s, 0.000158s

-- 57. Schur-elementary conjugate-Kostka route, cold/warm and reversed.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing QQ
time hallInnerProduct(S_{8,5,3}, e_{5,4,3,2,1,1})
time hallInnerProduct(S_{8,5,3}, e_{5,4,3,2,1,1})
time hallInnerProduct(e_{5,4,3,2,1,1}, S_{8,5,3})
-- 2026-07-11; codex QQ cold/warm/reversed: 0.000316s, 0.000411s, 0.000223s

-- 58. Registered diagonal pipeline on sparse and denser expansions.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time hallInnerProduct(Q_{8,4,2}, P_{8,4,2})
diagonalLeft = 2*Q_{8,4,2} + 3*Q_{7,5,2} - Q_{7,4,3} + 5*Q_{6,5,3}
diagonalRight = P_{8,4,2} - 2*P_{7,5,2} + 4*P_{7,4,3} + P_{6,5,3}
time hallInnerProduct(diagonalLeft, diagonalRight)
time hallInnerProduct(diagonalRight, diagonalLeft)
-- 2026-07-11; codex frac(QQ[t]) sparse/dense/reversed: 0.000201s, 0.000222s, 0.000212s

-- 59. Structured power-sum/Schur pipeline on sparse and denser expansions.
restart
needsPackage "SymmetricRings"
Sym = symmetricRing (frac(QQ[t]))
time hallInnerProduct(p_{8,4,2}, S_{8,4,2})
structuredPowerSums = 2*p_{8,4,2} + p_{7,5,2} - 3*p_{7,4,3} + p_{6,5,3}
structuredSchur = S_{8,4,2} - 2*S_{7,5,2} + S_{7,4,3} + 4*S_{6,5,3}
time hallInnerProduct(structuredPowerSums, structuredSchur)
time hallInnerProduct(structuredSchur, structuredPowerSums)
-- 2026-07-11; codex frac(QQ[t]) sparse/dense/reversed: 0.000319s, 0.006628s, 0.004429s
-- 2026-07-11 2:00pm EDT; standalone cold dense 3-run median: 0.004856s; runs: 0.005022s, 0.004856s, 0.004833s

-- 60. Unconditional fallback oracle. Run this block in a fresh process with:
-- M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_PIPELINE=fallback-power-sums M2
restart
needsPackage "SymmetricRings"
Sym = symmetricRing QQ
time hallInnerProduct(S_{8,5,3}, h_{5,4,3,2,1,1})
time hallInnerProduct(S_{8,5,3}, h_{5,4,3,2,1,1})
time hallInnerProduct(h_{5,4,3,2,1,1}, S_{8,5,3})
-- 2026-07-11; codex QQ forced fallback cold/warm/reversed: 0.004303s, 0.000836s, 0.000793s
