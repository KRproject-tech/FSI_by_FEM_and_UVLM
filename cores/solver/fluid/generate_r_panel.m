%% パネルノード点計算 [m]



%% 翼表面

%%[0-0] 点① [-]
r_panel_vec_1v = Sc_mat_panel_global_1*q_vec;
r_panel_vec_1 = reshape( r_panel_vec_1v, 3, []).';

%%[0-1] 点② [-]　(平板末端のパネルノード点は，平板末端ノード座標より補間)
r_panel_vec_2v = Sc_mat_panel_global_2*q_vec;
r_panel_vec_2 = reshape( r_panel_vec_2v, 3, []).';

%%[0-2] 点③ [-] (平板末端のパネルノード点は，平板末端ノード座標より補間)
r_panel_vec_3v = Sc_mat_panel_global_3*q_vec;
r_panel_vec_3 = reshape( r_panel_vec_3v, 3, []).';

%%[0-3] 点④ [-]
r_panel_vec_4v = Sc_mat_panel_global_4*q_vec;
r_panel_vec_4 = reshape( r_panel_vec_4v, 3, []).';



%%[0-4] 点① [-]
dt_r_panel_vec_1v = Sc_mat_panel_global_1*dt_q_vec;
dt_r_panel_vec_1 = reshape( dt_r_panel_vec_1v, 3, []).';

%%[0-5] 点② [-]　(平板末端のパネルノード点は，平板末端ノード座標より補間)
dt_r_panel_vec_2v = Sc_mat_panel_global_2*dt_q_vec;
dt_r_panel_vec_2 = reshape( dt_r_panel_vec_2v, 3, []).';

%%[0-6] 点③ [-] (平板末端のパネルノード点は，平板末端ノード座標より補間)
dt_r_panel_vec_3v = Sc_mat_panel_global_3*dt_q_vec;
dt_r_panel_vec_3 = reshape( dt_r_panel_vec_3v, 3, []).';

%%[0-7] 点④ [-]
dt_r_panel_vec_4v = Sc_mat_panel_global_4*dt_q_vec;
dt_r_panel_vec_4 = reshape( dt_r_panel_vec_4v, 3, []).';





%% 翼裏面

%%[0-0] 点① [-]
r_panel_vec_1v_1 = Sc_mat_panel_global_1*q_vec_1;
r_panel_vec_1_1 = reshape( r_panel_vec_1v_1, 3, []).';

%%[0-1] 点② [-]　(平板末端のパネルノード点は，平板末端ノード座標より補間)
r_panel_vec_2v_1 = Sc_mat_panel_global_2*q_vec_1;
r_panel_vec_2_1 = reshape( r_panel_vec_2v_1, 3, []).';

%%[0-2] 点③ [-] (平板末端のパネルノード点は，平板末端ノード座標より補間)
r_panel_vec_3v_1 = Sc_mat_panel_global_3*q_vec_1;
r_panel_vec_3_1 = reshape( r_panel_vec_3v_1, 3, []).';

%%[0-3] 点④ [-]
r_panel_vec_4v_1 = Sc_mat_panel_global_4*q_vec_1;
r_panel_vec_4_1 = reshape( r_panel_vec_4v_1, 3, []).';



%%[0-4] 点① [-]
dt_r_panel_vec_1v_1 = Sc_mat_panel_global_1*dt_q_vec_1;
dt_r_panel_vec_1_1 = reshape( dt_r_panel_vec_1v_1, 3, []).';

%%[0-5] 点② [-]　(平板末端のパネルノード点は，平板末端ノード座標より補間)
dt_r_panel_vec_2v_1 = Sc_mat_panel_global_2*dt_q_vec_1;
dt_r_panel_vec_2_1 = reshape( dt_r_panel_vec_2v_1, 3, []).';

%%[0-6] 点③ [-] (平板末端のパネルノード点は，平板末端ノード座標より補間)
dt_r_panel_vec_3v_1 = Sc_mat_panel_global_3*dt_q_vec_1;
dt_r_panel_vec_3_1 = reshape( dt_r_panel_vec_3v_1, 3, []).';

%%[0-7] 点④ [-]
dt_r_panel_vec_4v_1 = Sc_mat_panel_global_4*dt_q_vec_1;
dt_r_panel_vec_4_1 = reshape( dt_r_panel_vec_4v_1, 3, []).';




%% コロケーション点計算 [m]
%%% 1:要素番号，2:座標成分    (翼表面)
rc_vec = Sc_mat_col_global*q_vec;
rc_vec = reshape( rc_vec, 3, []).';


%%% 1:要素番号，2:座標成分    (翼裏面)
rc_vec_1 = Sc_mat_col_global*q_vec_1;
rc_vec_1 = reshape( rc_vec_1, 3, []).';


%% コロケーション点の変位速度 [-]


%%% 1:要素番号，2:座標成分    (翼表面)
dt_rc_vec = Sc_mat_col_global*dt_q_vec;
dt_rc_vec = reshape( dt_rc_vec, 3, []).';

%%% 1:要素番号，2:座標成分    (翼裏面)
dt_rc_vec_1 = Sc_mat_col_global*dt_q_vec_1;
dt_rc_vec_1 = reshape( dt_rc_vec_1, 3, []).';




%% 統合

r_panel_vec_1_all = [   r_panel_vec_1;
                        r_panel_vec_1_1];
r_panel_vec_2_all = [   r_panel_vec_2;
                        r_panel_vec_2_1];
r_panel_vec_3_all = [   r_panel_vec_3;
                        r_panel_vec_3_1];
r_panel_vec_4_all = [   r_panel_vec_4;
                        r_panel_vec_4_1];
                    
dt_r_panel_vec_1_all = [    dt_r_panel_vec_1;
                            dt_r_panel_vec_1_1];
dt_r_panel_vec_2_all = [    dt_r_panel_vec_2;
                            dt_r_panel_vec_2_1];
dt_r_panel_vec_3_all = [    dt_r_panel_vec_3;
                            dt_r_panel_vec_3_1];
dt_r_panel_vec_4_all = [    dt_r_panel_vec_4;
                            dt_r_panel_vec_4_1];
                    

rc_vec_all = [  rc_vec;
                rc_vec_1];

dt_rc_vec_all = [   dt_rc_vec;
                    dt_rc_vec_1];
