(11:48 2016/02/01)
[1] 基本セット (GUI.figを起動する)

(10:33 2017/06/19)
[1] 構造減衰があるときは，膜ひずみ剛性において，Qem = Qem(q_k) + 1/2*Δt*dq_Qem(q_k)dt_qを用いる．(エネルギ収支の精度向上)

(13:16 2017/11/16)
[1] 法線ベクトルの勾配∂q_(n/||n||)の計算方法を修正．
[2] wakeモデルの時間発展にRK4法を使用．流体力計算で反復計算を行う．

(23:33 2019/06/03)
[1] 計算高速化・計算時間の算出

(12:07 2023/11/20)
[1] Plate_FEM_VLM_NONDIM_10 をベースに，弱連成・強連成 解法の選択を追加．
※人工付加質量法[1][2]では，人工付加質量行列の成分が実際の付加質量に近くないと収束が難しいため，自由度の大きいFEMで人工付加質量行列を決定するのは困難である．
[1] F. Be´langer, M.P. Pa?¨doussis, E. de Langre, Time-marching analysis of fluid-coupled system with large added mass, AIAA Journal, Vol. 33, (1995).
[2] C. Yvin et al., Added mass evaluation with a finite-volume solver for applications in fluid?structure interaction problems solved with co-simulation, Journal of Fluids and Structures, Vol. 81, (2018) .

(16:10 2024/03/21)
[1] 平板スパン方向中央の，自由端変位の時刻歴・変位のスナップショット，流速分布の等高線のplotを追加．