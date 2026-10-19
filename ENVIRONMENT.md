# 环境与复现

## Python 后处理

Python 3.11 或更新版本，依赖见 [requirements.txt](requirements.txt)：NumPy、Matplotlib 和 Pillow。发布前在 Python 3.11 下执行文件与科学数值核验、表格导出、八张图和两段动画重建。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
python reproduce.py --figures --animations
python -m unittest discover -s scripts -p 'test_*.py'
```

默认命令从现有文本输出、fchk 和已审阅数据离线工作，不调用电子结构软件。`--output <新目录>` 可指定独立输出位置，已有目录会被拒绝，以免覆盖文件。

导出表格：四方法电子能比较、16 点 PM7 扫描、驻点热化学、56 点 IRC 自旋与键长。核验同时检查原始失败扫描状态、源/公开文件校验值、Hessian 投影的 27 个振动、稳定性、同几何单点、CC 收敛/诊断及 Mulliken/Hirshfeld 全原子表。

## 报告排版

使用 XeLaTeX，需 `ctex`、Fandol 字体、TikZ、booktabs、tabularx 等 TeX 源中列明的包；完整 TeX Live/MacTeX 可提供。

```bash
python reproduce.py --tex
```

该命令连续编译两遍，在本次独立输出目录生成 PDF 和日志，不覆盖已发布 PDF。报告日期固定在 TeX 中，不依赖重新编译日期。

## 原计算软件

主计算使用 Gaussian 16 Rev. A.03；耦合簇对照使用 ORCA 6.1.1。仓库仅提供本次输入/文本输出/格式化检查点，不分发这些程序，也不自动配置或启动计算。沿 checkpoint 衔接的 Gaussian 输入保留历史科学设置；若重新计算，需依据前序任务重新生成所需 checkpoint，不能把后处理脚本当作提交工具。

## 复现判据

能量、频率、电子数、自旋布居及源文件校验值为主要核验对象。重绘图像可能随 Matplotlib、字体或平台产生排版差异，不要求重建图片/PDF 二进制完全一致。发布版 PDF 的校验值保存在 `analysis/report-review-20261019/final-report-audit.json`，其 16 页已逐页检查。
