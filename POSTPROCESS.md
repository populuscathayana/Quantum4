# 后处理范围

`reproduce.py` 先核验公开副本，再向新目录导出数据；`--figures` 重建 PM7 完整曲线、IRC 能量/自旋、IRC 键长/⟨S²⟩、电子能/Gibbs 阶梯图、热化学分解、四方法比较、TS 模式和 IRC 关键结构，共八张图。`--animations` 还重建 32 帧 TS 模式与 56 帧 IRC 结构动画。

## 能量与路径

- PM7 与 DFT 单独取零点，不混成同一能量曲线。
- IRC 采用 Gaussian 接受的校正后几何，不用预测点代替；路径坐标单位为 √amu·bohr。
- 离散 E/G 阶梯的虚线只是连接驻点的示意线，不是连续自由能路径。
- 四种方法的固定几何比较分别以各自 R 的电子能为零。ORCA 与 Gaussian 的绝对能量不相减。
- 298.15 K、1 atm 的 RRHO 热化学由原 DFT R/TS/P 频率计算给出；三者转动对称数均为 1。

## 结构、模式与密度

TS 位移来自公开 fchk 的 Cartesian Hessian，经质量加权并投影掉六个刚体自由度；动画振幅只用于展示。IRC 动画使用真实接受点，保持原子身份，仅作平移和正旋转配准。

已发布的 `experiment-report/figures/spin-density.png` 是真实 SCF α−β 密度的 ±0.02 e a₀⁻³ 等值面，不是轨道示意图。网格、积分、Spin 与 FAlpha−FBeta 逐点差异核验保存在 `analysis/spin-density-20261019/audit.json`。

为控制仓库体积，Spin/FAlpha/FBeta cube 网格不上传；默认离线复现保留这张已核验成图，不重新生成密度网格。R、TS、P 的完整 fchk 均已公开，分别位于 `data/calculations/reactant-opt/`、`ts-opt/`、`product-endpoint-opt/`。如需重建网格，可使用合法安装的 Gaussian `cubegen` 按核验文件中的同一原点、轴向步长和网格数导出 `Spin=SCF`、`FAlpha=SCF`、`FBeta=SCF`，先验证积分及逐点 α−β 一致，再绘制等值面。

## 数据映射

`data/manifest.json` 的 `source` 是原项目相对路径，`path` 是公开副本路径；`original_sha256` 对应历史科学核验，`sha256` 对应公开文件。路径脱敏不改变科学数据。未分发的原始 cube、检查点、调度/回收记录仍可在历史核验中被引用，但默认后处理不依赖这些文件。
