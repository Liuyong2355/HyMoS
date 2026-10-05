function schema = hymos_form_schema()
% UI metadata only: actual parsing and solving remain in the original chain.
schema.groups = {'物理参数', '网格与阶数', '求解方法', '多重网格', '完整配置'};
% group, INI key, Chinese label, explanation, dropdown items (empty = text)
schema.rows = {
1, 'Const/Kn', 'Knudsen 数', '稀薄度（>0）。', {};
1, 'LeftWall/Temperature', '左壁温度 T_L', '正数。', {};
1, 'RightWall/Temperature', '右壁温度 T_R', '正数。', {};
1, 'LeftWall/TangentialVelocity', '左壁速度 u2_L', '有限实数；仅切向。', {};
1, 'RightWall/TangentialVelocity', '右壁速度 u2_R', '有限实数；仅切向。', {};
2, 'Mesh/Nx', '网格单元数', '整数，默认 128。', {};
2, 'Mesh/L0', '区间左端点', '默认 -0.5。', {};
2, 'Mesh/Len', '区间长度', '正数，默认 1。', {};
2, 'System/ORDER', '矩展开阶数', '越高越耗时。', {};
3, 'Solver/SolverType', '求解框架', 'single 或 nmg。', {'single', 'nmg'};
3, 'Solver/FIM', 'FIM 方法', '当前默认 3。', {'1', '2', '3'};
3, 'Error/Tol', '残差容差', '默认 1e-8。', {};
3, 'Time/CFL', '矩方程 CFL', '默认 0.8。', {};
3, 'Time/CFL_HEs', '宏观方程 CFL', '默认 0.8。', {};
3, 'SingleGridMethod/EulerEqsSolvingSteps', 'Euler 迭代步数', '默认 600。', {};
3, 'System/n_thread', '计算线程数', '按本机资源设置。', {};
4, 'MultiGridMethod/PreSmoothingSteps', '前平滑次数', '默认 2。', {};
4, 'MultiGridMethod/PostSmoothingSteps', '后平滑次数', '默认 2。', {};
4, 'MultiGridMethod/CoarsestGridSmoothingSteps', '最粗层平滑次数', '默认 4。', {};
4, 'MultiGridMethod/CoarsestGridSize', '最粗网格单元数', '默认 8。', {};
4, 'MultiGridMethod/CycleType', '网格循环', '1：V；2：W；3：F。', {'1', '2', '3'}
};
end
