# Edge Computing Codes
- [Edge Computing Codes](#Edge-Computing-Codes)
  - [UIC18](#UIC18)
    - [LODCO algorithm](#LODCO-algorithm)
    - [LODCO-based Greedy algorithms](#LODCO-based-Greedy-algorithms)
  - [Final Notes](#Final-Notes)
    - [About the Author](#About-the-Author)

Edge computing is the practice of processing data near the edge of the network, where the data is being generated, instead of in a centralised data-processing warehouse. This repos provides a better implementation of proposed algorithms in Edge Computing (it may not the same as the original algorithm).

## UIC18
The folder `UIC18` includes codes for the following paper:

**Hailiang Zhao**, Wei Du, Wei Liu, Tao Lei and Qiwang Lei, *QoE Aware and Cell Capacity Enhanced Computation Offloading for Multi-Server Mobile Edge Computing Systems with Energy Harvesting Devices.* In: **Proceedings of the 15th IEEE International Conference on Ubiquitous Intelligence and Computing (UIC'18)**, Guangzhou, China, 2018.

### LODCO algorithm
This paper is based on **the Lyapunov Optimization-based Dynamic Computation Offloading (LODCO) algorithm**, which was proposed in the following paper:

Y. Mao, J. Zhang and K. B. Letaief, *Dynamic Computation Offloading for Mobile-Edge Computing With Energy Harvesting Devices*. In: **IEEE Journal on Selected Areas in Communications**, vol. 34, no. 12, pp. 3590-3605, Dec. 2016.

We implemented **the LODCO algorithm** (`LODCO.m`), and the simulation results are as follows (50000 time slots).
<div align=center>
    <img src="./UIC18/figures/LODCO_battery.svg" width="450"/><img src="./UIC18/figures/LODCO_cost.svg" width="450"/>
</div>
<div align=center>
<img src="./UIC18/figures/LODCO_ratio.svg" width="800"/>
</div>

### LODCO-based Greedy algorithms
**LODCO-based Greedy algorithm** was proposed for computation offloading in *multi-user multi-server* Mobile Edge Computing (MEC) systems. We implemented the algorithm (`greedy_LODCO.m`). The simulation results are as follows (2000 time slots, 10 mobile devices and 8 MEC servers).
<div align=center>
    <img src="./UIC18/figures/10devices_battery.svg" width="450"/><img src="./UIC18/figures/10devices_cost.svg" width="450"/>
</div>
<div align=center>
    <img src="./UIC18/figures/10devices_ratio.svg" width="800"/>
</div>

**LODCO-based epsilon-Greedy algorithm** was proposed to increase the ratio of MEC server execution (`eps_greedy_LODCO.m`). The ratio comparasion is as follows. In the left figure (above), epsilon is set as 0.3, and in the right one (below), epsilon is 0.8. The ratio of MEC server execution increases as epsilon increases.
<div align=center>
<img src="./UIC18/figures/10devices_ratio_eps03.svg" width="450"/><img src="./UIC18/figures/10devices_ratio_eps08.svg" width="450"/>
</div>

**LODCO-based Genetic Algorithm with Greedy Policy** was an alternative way for computation offloading in *multi-user multi-server** MEC systems. The code will be added later.

## Final Notes
If you have used the codes in your research works, we would appreciate citation to the paper mentioned before:

**Hailiang Zhao**, Wei Du, Wei Liu, Tao Lei and Qiwang Lei, *QoE Aware and Cell Capacity Enhanced Computation Offloading for Multi-Server Mobile Edge Computing Systems with Energy Harvesting Devices.* In: **Proceedings of the 15th IEEE International Conference on Ubiquitous Intelligence and Computing (UIC'18)**, Guangzhou, China, 2018.

### About the Author
**Hailiang Zhao** @ [CCNT Lab](http://www.cs.zju.edu.cn/kejizhan/lab/CCNT.html), [ZJU-CS](http://www.cs.zju.edu.cn/)
* Email: hliangzhao97 {AT} gmail.com
* GitHub: https://github.com/hliangzhao
* Profile: https://hliangzhao.github.io/CV/