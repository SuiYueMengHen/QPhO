# 从电磁场拉格朗日量推导麦克斯韦方程组

## 第一篇 必备数学基础
### 第一章 矢量代数基础
#### 1.1 基本定义
**定义1.1.1 三维矢量**：在三维欧几里得空间中，矢量$\boldsymbol{F}$可表示为直角坐标系下的分量形式：
$$\boldsymbol{F} = F_x \boldsymbol{e}_x + F_y \boldsymbol{e}_y + F_z \boldsymbol{e}_z$$
其中$\boldsymbol{e}_x,\boldsymbol{e}_y,\boldsymbol{e}_z$为x、y、z轴方向的单位正交基矢量，满足$\boldsymbol{e}_i \cdot \boldsymbol{e}_j = \begin{cases}1, & i=j \\ 0, & i\neq j\end{cases}$。

为简化书写，约定：拉丁字母指标$i,j,k,l,m,n$仅取1,2,3，分别对应x,y,z轴，即$\boldsymbol{e}_1=\boldsymbol{e}_x,\boldsymbol{e}_2=\boldsymbol{e}_y,\boldsymbol{e}_3=\boldsymbol{e}_z$，矢量可简写为$\boldsymbol{F} = \sum_{i=1}^3 F_i \boldsymbol{e}_i$。

#### 1.2 克罗内克δ函数
**定义1.2.1 克罗内克δ函数（Kronecker Delta）**：
$$\delta_{ij} = \begin{cases}1, & i=j \\ 0, & i\neq j\end{cases}, \quad i,j=1,2,3$$

**性质1.2.1 筛选性**：对任意序列$a_j$，有$\sum_{j=1}^3 \delta_{ij} a_j = a_i$
**证明**：当$j=i$时，$\delta_{ii}=1$，项为$a_i$；当$j\neq i$时，$\delta_{ij}=0$，项为0，求和后仅余$a_i$，证毕。

#### 1.3 列维-奇维塔符号
**定义1.3.1 列维-奇维塔符号（Levi-Civita Symbol）**：三维全反对称符号，定义为
$$\varepsilon_{ijk} = \begin{cases}1, & (i,j,k)为(1,2,3)的偶置换 \\ -1, & (i,j,k)为(1,2,3)的奇置换 \\ 0, & 存在任意两个指标相等\end{cases}$$
例：$\varepsilon_{123}=1$，$\varepsilon_{213}=-1$（交换1和2，奇置换），$\varepsilon_{112}=0$（两个指标重复）。

**性质1.3.1 全反对称性**：交换任意两个指标，符号取反，即$\varepsilon_{ijk}=-\varepsilon_{jik}=-\varepsilon_{ikj}=-\varepsilon_{kji}$

**性质1.3.2 核心恒等式（ε-δ恒等式）**：
$$\sum_{i=1}^3 \varepsilon_{ijk}\varepsilon_{imn} = \delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km}$$
**证明**：取$j=2,k=3,m=2,n=3$，左边$=\varepsilon_{123}\varepsilon_{123}+\varepsilon_{223}\varepsilon_{223}+\varepsilon_{323}\varepsilon_{323}=1\times1+0+0=1$，右边$=\delta_{22}\delta_{33}-\delta_{23}\delta_{32}=1\times1-0=1$，成立；其余取值可逐一验证，证毕。

#### 1.4 矢量的点积与叉积
**定义1.4.1 点积（标量积）**：两个矢量$\boldsymbol{A},\boldsymbol{B}$的点积定义为
$$\boldsymbol{A}\cdot\boldsymbol{B} = |\boldsymbol{A}||\boldsymbol{B}|\cos\theta$$
其中$\theta$为两矢量的夹角。

**定理1.4.1 点积的分量形式**：$\boldsymbol{A}\cdot\boldsymbol{B} = \sum_{i=1}^3 A_i B_i$
**证明**：$\boldsymbol{A}\cdot\boldsymbol{B} = \left(\sum_{i=1}^3 A_i \boldsymbol{e}_i\right)\cdot\left(\sum_{j=1}^3 B_j \boldsymbol{e}_j\right) = \sum_{i,j=1}^3 A_i B_j (\boldsymbol{e}_i\cdot\boldsymbol{e}_j) = \sum_{i,j=1}^3 A_i B_j \delta_{ij} = \sum_{i=1}^3 A_i B_i$，证毕。

**定义1.4.2 叉积（矢量积）**：两个矢量$\boldsymbol{A},\boldsymbol{B}$的叉积定义为
$$\boldsymbol{A}\times\boldsymbol{B} = |\boldsymbol{A}||\boldsymbol{B}|\sin\theta \cdot \boldsymbol{e}_n$$
其中$\boldsymbol{e}_n$为垂直于$\boldsymbol{A},\boldsymbol{B}$所张平面的单位矢量，方向满足右手螺旋定则。

**定理1.4.2 叉积的分量形式**：$(\boldsymbol{A}\times\boldsymbol{B})_k = \sum_{i,j=1}^3 \varepsilon_{ijk} A_i B_j$
**证明**：以$k=1$为例，$(\boldsymbol{A}\times\boldsymbol{B})_1 = A_2 B_3 - A_3 B_2$，右边$\sum_{i,j=1}^3 \varepsilon_{ij1}A_i B_j = \varepsilon_{231}A_2 B_3 + \varepsilon_{321}A_3 B_2 = 1\cdot A_2 B_3 + (-1)\cdot A_3 B_2$，与左边相等；$k=2,k=3$可同理验证，证毕。

---

### 第二章 矢量分析与场论基础
#### 2.1 场的定义
**定义2.1.1 标量场**：若三维空间中每一点$(x,y,z)$和时刻$t$，都对应一个标量$\phi(x,y,z,t)$，则称$\phi$为标量场，例如温度场、电势场。

**定义2.1.2 矢量场**：若三维空间中每一点$(x,y,z)$和时刻$t$，都对应一个矢量$\boldsymbol{F}(x,y,z,t)$，则称$\boldsymbol{F}$为矢量场，例如电场、磁场、流速场。

#### 2.2 梯度、散度、旋度的定义与分量形式
**定义2.2.1 哈密顿算子（Nabla算子）**：定义为兼具矢量性和微分性的算子
$$\nabla = \boldsymbol{e}_1 \frac{\partial}{\partial x_1} + \boldsymbol{e}_2 \frac{\partial}{\partial x_2} + \boldsymbol{e}_3 \frac{\partial}{\partial x_3} = \sum_{i=1}^3 \boldsymbol{e}_i \frac{\partial}{\partial x_i}$$

**定义2.2.2 梯度**：对标量场$\phi$，梯度定义为
$$\nabla\phi = \sum_{i=1}^3 \boldsymbol{e}_i \frac{\partial \phi}{\partial x_i}$$
梯度是矢量，方向为标量场变化最快的方向，大小为最大变化率。

**定义2.2.3 散度**：对矢量场$\boldsymbol{F}$，散度定义为
$$\nabla\cdot\boldsymbol{F} = \sum_{i=1}^3 \frac{\partial F_i}{\partial x_i}$$
散度是标量，描述矢量场在某点的「源」或「汇」的强度。

**定义2.2.4 旋度**：对矢量场$\boldsymbol{F}$，旋度定义为
$$(\nabla\times\boldsymbol{F})_k = \sum_{i,j=1}^3 \varepsilon_{ijk} \frac{\partial F_j}{\partial x_i}$$
旋度是矢量，描述矢量场在某点的旋转特性。

#### 2.3 核心积分定理（高数内容回顾与推广）
**定理2.3.1 高斯散度定理（高数高斯公式）**
设$V$为三维空间中的有界闭区域，其边界$S$为分片光滑的闭合曲面，矢量场$\boldsymbol{F}$在$V$上具有连续的一阶偏导数，则
$$\int_V \nabla\cdot\boldsymbol{F} \, dV = \oint_S \boldsymbol{F}\cdot d\boldsymbol{S}$$
其中$d\boldsymbol{S} = \boldsymbol{e}_n dS$，$\boldsymbol{e}_n$为曲面$S$的外法线单位矢量，$dS$为面元面积。
物理意义：矢量场在闭合曲面上的通量，等于场在曲面所围体积内的散度的体积分。

**定理2.3.2 斯托克斯旋度定理（高数斯托克斯公式）**
设$\Sigma$为三维空间中分片光滑的有向曲面，其边界$L$为分段光滑的闭合有向曲线，$L$的正向与$\Sigma$的法向满足右手螺旋定则，矢量场$\boldsymbol{F}$在$\Sigma$上具有连续的一阶偏导数，则
$$\int_\Sigma (\nabla\times\boldsymbol{F})\cdot d\boldsymbol{S} = \oint_L \boldsymbol{F}\cdot d\boldsymbol{l}$$
其中$d\boldsymbol{l}$为曲线$L$的线元矢量，方向沿$L$的正向。
物理意义：矢量场沿闭合曲线的环量，等于场在曲线所围曲面上的旋度的通量。

#### 2.4 常用矢量恒等式（全证明）
**恒等式2.4.1** 对任意标量场$\phi$，有$\nabla\times(\nabla\phi) = 0$
**证明**：$(\nabla\times(\nabla\phi))_k = \sum_{i,j=1}^3 \varepsilon_{ijk} \frac{\partial}{\partial x_i}\left( \frac{\partial \phi}{\partial x_j} \right) = \sum_{i,j=1}^3 \varepsilon_{ijk} \frac{\partial^2 \phi}{\partial x_i \partial x_j}$
由高数混合偏导数相等定理，$\frac{\partial^2 \phi}{\partial x_i \partial x_j} = \frac{\partial^2 \phi}{\partial x_j \partial x_i}$，交换指标$i,j$得：
$$(\nabla\times(\nabla\phi))_k = \sum_{i,j=1}^3 \varepsilon_{jik} \frac{\partial^2 \phi}{\partial x_j \partial x_i} = -\sum_{i,j=1}^3 \varepsilon_{ijk} \frac{\partial^2 \phi}{\partial x_i \partial x_j} = -(\nabla\times(\nabla\phi))_k$$
因此$2(\nabla\times(\nabla\phi))_k=0$，即$\nabla\times(\nabla\phi)=0$，证毕。

**恒等式2.4.2** 对任意矢量场$\boldsymbol{F}$，有$\nabla\cdot(\nabla\times\boldsymbol{F}) = 0$
**证明**：$\nabla\cdot(\nabla\times\boldsymbol{F}) = \sum_{k=1}^3 \frac{\partial}{\partial x_k} (\nabla\times\boldsymbol{F})_k = \sum_{k,i,j=1}^3 \varepsilon_{ijk} \frac{\partial^2 F_j}{\partial x_k \partial x_i}$
同理，交换$k,i$，利用混合偏导数相等和ε的反对称性，可得该式等于自身的相反数，因此为0，证毕。

---

### 第三章 变分法基础（核心前置知识）
#### 3.1 泛函与变分的定义
**定义3.1.1 泛函**：设$C$为满足一定条件的函数集合，若对集合中任意一个函数$y(x)$，都有一个实数$J$与之对应，则称$J$为函数$y(x)$的泛函，记为$J[y(x)]$。
通俗理解：泛函是「函数的函数」，自变量是函数，因变量是实数，例如定积分$J[y] = \int_a^b y(x) dx$就是一个泛函。

**定义3.1.2 函数的变分**：对函数$y(x)$，其微小的、虚的变动$\delta y(x) = \tilde{y}(x) - y(x)$称为函数$y(x)$的变分，其中$\tilde{y}(x)$是与$y(x)$非常接近的函数。

变分的核心性质（与微分完全一致）：
1.  变分与微分可交换顺序：$\delta \left( \frac{dy}{dx} \right) = \frac{d}{dx}(\delta y)$
2.  变分与积分可交换顺序：$\delta \int_a^b y(x) dx = \int_a^b \delta y(x) dx$
3.  线性运算法则：$\delta(y_1+y_2)=\delta y_1+\delta y_2$，$\delta(y_1 y_2)=y_1 \delta y_2 + y_2 \delta y_1$

**定义3.1.3 泛函的变分**：泛函$J[y]$的增量$\Delta J = J[y+\delta y] - J[y]$可展开为
$$\Delta J = L[y,\delta y] + o(\delta y)$$
其中$L[y,\delta y]$是关于$\delta y$的线性项，$o(\delta y)$是$\delta y$的高阶无穷小项，则称线性项$L[y,\delta y]$为泛函$J[y]$的变分，记为$\delta J$。

#### 3.2 最小作用量原理与一元函数泛函的欧拉-拉格朗日方程
**定理3.2.1 最小作用量原理（哈密顿原理）**：力学系统的真实运动，对应作用量泛函$S$取极值，即$\delta S = 0$，且在时间的端点$t_1,t_2$处，广义坐标的变分为0：$\delta q(t_1)=\delta q(t_2)=0$。

对于单个质点的一维运动，作用量泛函定义为
$$S[q(t)] = \int_{t_1}^{t_2} L(q(t), \dot{q}(t), t) dt$$
其中$q(t)$为广义坐标，$\dot{q}(t)=\frac{dq}{dt}$为广义速度，$L$为拉格朗日量，是$q,\dot{q},t$的函数。

**定理3.2.2 一元函数泛函的欧拉-拉格朗日方程**
若泛函$S[q(t)] = \int_{t_1}^{t_2} L(q,\dot{q},t) dt$取极值，且端点变分$\delta q(t_1)=\delta q(t_2)=0$，则极值函数$q(t)$满足
$$\frac{d}{dt}\left( \frac{\partial L}{\partial \dot{q}} \right) - \frac{\partial L}{\partial q} = 0$$
**证明**：
对作用量取变分，交换变分与积分顺序：
$$\delta S = \delta \int_{t_1}^{t_2} L dt = \int_{t_1}^{t_2} \delta L dt$$
由多元函数的全微分法则，$L$的变分为：
$$\delta L = \frac{\partial L}{\partial q} \delta q + \frac{\partial L}{\partial \dot{q}} \delta \dot{q}$$
代入变分表达式：
$$\delta S = \int_{t_1}^{t_2} \left( \frac{\partial L}{\partial q} \delta q + \frac{\partial L}{\partial \dot{q}} \delta \dot{q} \right) dt$$
由变分的性质$\delta \dot{q} = \frac{d}{dt}(\delta q)$，对第二项做分部积分（高数分部积分法：$\int_a^b u dv = uv|_a^b - \int_a^b v du$，令$u=\frac{\partial L}{\partial \dot{q}}$，$dv=\frac{d}{dt}(\delta q) dt$，则$du=\frac{d}{dt}(\frac{\partial L}{\partial \dot{q}})dt$，$v=\delta q$）：
$$\int_{t_1}^{t_2} \frac{\partial L}{\partial \dot{q}} \frac{d}{dt}(\delta q) dt = \left. \frac{\partial L}{\partial \dot{q}} \delta q \right|_{t_1}^{t_2} - \int_{t_1}^{t_2} \frac{d}{dt}\left( \frac{\partial L}{\partial \dot{q}} \right) \delta q dt$$
由端点条件$\delta q(t_1)=\delta q(t_2)=0$，第一项边界项为0，因此：
$$\delta S = \int_{t_1}^{t_2} \left[ \frac{\partial L}{\partial q} - \frac{d}{dt}\left( \frac{\partial L}{\partial \dot{q}} \right) \right] \delta q dt$$
泛函取极值的条件是$\delta S=0$，且$\delta q$是任意的（除端点外），因此被积函数必须恒为0，移项后即得欧拉-拉格朗日方程，证毕。

#### 3.3 场的欧拉-拉格朗日方程（推广到多元函数）
对于连续场，场量$\phi(x,y,z,t)$是时空坐标的多元函数，相当于无穷多个广义坐标（每个空间点的场量都是一个广义坐标）。此时，拉格朗日量$L$是拉格朗日密度$\mathcal{L}$对全空间的体积分：
$$L(t) = \int_{V} \mathcal{L}(\phi, \partial_t \phi, \partial_x \phi, \partial_y \phi, \partial_z \phi, x,y,z,t) dV$$
其中$\mathcal{L}$是场量$\phi$、场量的所有一阶偏导数、时空坐标的函数。

作用量泛函为拉格朗日量对时间的积分，即拉格朗日密度对四维时空的积分：
$$S[\phi] = \int_{t_1}^{t_2} L(t) dt = \int_{t_1}^{t_2} \int_V \mathcal{L} dV dt = \int_\Omega \mathcal{L} d^4x$$
其中$\Omega$为四维时空积分区域，$d^4x = dt dx dy dz$为四维时空体积元。

**定义3.3.1 爱因斯坦求和约定**（后续全程使用）：
1.  希腊字母指标$\mu,\nu,\lambda,\sigma$取0,1,2,3，对应四维时空坐标$x^0=ct, x^1=x, x^2=y, x^3=z$；
2.  当一个指标在一项中同时出现一次上标、一次下标时，自动对该指标的所有取值求和，无需书写求和号；
3.  重复出现的求和指标称为「哑指标」，可任意更换字母；仅出现一次的指标称为「自由指标」，等式两边的自由指标必须完全一致。

例：$\partial_\mu \left( \frac{\partial \mathcal{L}}{\partial (\partial_\mu \phi)} \right) = \sum_{\mu=0}^3 \partial_\mu \left( \frac{\partial \mathcal{L}}{\partial (\partial_\mu \phi)} \right)$，其中$\partial_0 = \frac{\partial}{\partial x^0} = \frac{1}{c}\frac{\partial}{\partial t}$，$\partial_1=\frac{\partial}{\partial x},\partial_2=\frac{\partial}{\partial y},\partial_3=\frac{\partial}{\partial z}$。

**定理3.3.1 场的欧拉-拉格朗日方程**
若场的作用量$S[\phi]$取极值，且在时空区域$\Omega$的边界上，场量的变分$\delta\phi|_{\text{边界}}=0$，则场量$\phi$满足
$$\partial_\mu \left( \frac{\partial \mathcal{L}}{\partial (\partial_\mu \phi)} \right) - \frac{\partial \mathcal{L}}{\partial \phi} = 0$$
**证明**：
对作用量取变分，交换变分与积分顺序：
$$\delta S = \int_\Omega \delta \mathcal{L} d^4x$$
$\mathcal{L}$的变分为全微分：
$$\delta \mathcal{L} = \frac{\partial \mathcal{L}}{\partial \phi} \delta\phi + \frac{\partial \mathcal{L}}{\partial (\partial_t \phi)} \delta(\partial_t \phi) + \frac{\partial \mathcal{L}}{\partial (\partial_x \phi)} \delta(\partial_x \phi) + \frac{\partial \mathcal{L}}{\partial (\partial_y \phi)} \delta(\partial_y \phi) + \frac{\partial \mathcal{L}}{\partial (\partial_z \phi)} \delta(\partial_z \phi)$$
由变分与偏导数可交换顺序：$\delta(\partial_\mu \phi) = \partial_\mu (\delta\phi)$，对每个偏导数项做分部积分，利用边界条件$\delta\phi|_{\text{边界}}=0$，所有边界项均为0，最终得到：
$$\delta S = \int_\Omega \left[ \frac{\partial \mathcal{L}}{\partial \phi} - \partial_\mu\left( \frac{\partial \mathcal{L}}{\partial (\partial_\mu \phi)} \right) \right] \delta\phi d^4x$$
由$\delta S=0$且$\delta\phi$任意，被积函数必须恒为0，移项后即得场的欧拉-拉格朗日方程，证毕。

---

## 第二篇 电磁场的基本物理量与拉格朗日密度
### 第四章 电磁场的基本物理量定义
#### 4.1 电荷与电流
**定义4.1.1 电荷密度$\rho(x,y,z,t)$**：单位体积内的电荷量，单位为$\text{C/m}^3$。
**定义4.1.2 电流密度$\boldsymbol{J}(x,y,z,t)$**：单位时间内通过垂直于电流方向的单位面积的电荷量，单位为$\text{A/m}^2$，方向为正电荷的运动方向。

**定理4.1.1 电荷守恒定律（电流连续性方程）**
电荷既不能被创造，也不能被消灭，微分形式为：
$$\frac{\partial \rho}{\partial t} + \nabla\cdot\boldsymbol{J} = 0$$
**证明**：由电荷守恒，闭合曲面$S$内的总电荷量$Q=\int_V \rho dV$的减少率等于流出$S$的电流$I=\oint_S \boldsymbol{J}\cdot d\boldsymbol{S}$，即
$$-\frac{dQ}{dt} = \oint_S \boldsymbol{J}\cdot d\boldsymbol{S}$$
左边交换求导与积分顺序：$-\frac{dQ}{dt} = -\int_V \frac{\partial \rho}{\partial t} dV$；右边由高斯散度定理：$\oint_S \boldsymbol{J}\cdot d\boldsymbol{S} = \int_V \nabla\cdot\boldsymbol{J} dV$，因此
$$-\int_V \frac{\partial \rho}{\partial t} dV = \int_V \nabla\cdot\boldsymbol{J} dV$$
移项得$\int_V \left( \frac{\partial \rho}{\partial t} + \nabla\cdot\boldsymbol{J} \right) dV = 0$，由于体积$V$是任意的，被积函数必须恒为0，证毕。

#### 4.2 电磁场的标势与矢势
实验表明，电场$\boldsymbol{E}$和磁场$\boldsymbol{B}$可由标量势（标势）$\phi$和矢量势（矢势）$\boldsymbol{A}$完全描述，二者的定义如下：

**定义4.2.1 标势$\phi(x,y,z,t)$**：单位正电荷的电势能，单位为$\text{V}$（伏特）。
**定义4.2.2 矢势$\boldsymbol{A}(x,y,z,t)$**：描述磁场的矢量势，单位为$T·m$（特斯拉·米）。

**定理4.2.1 电场、磁场与势的基本关系**
$$\boldsymbol{E} = -\nabla\phi - \frac{\partial \boldsymbol{A}}{\partial t} \tag{4.2.1}$$
$$\boldsymbol{B} = \nabla\times\boldsymbol{A} \tag{4.2.2}$$
**证明**：
1.  对(4.2.2)式，两边取散度，由恒等式2.4.2，$\nabla\cdot(\nabla\times\boldsymbol{A})=0$，即$\nabla\cdot\boldsymbol{B}=0$，与磁场无磁单极子的实验规律一致；
2.  对(4.2.1)式，两边取旋度，由恒等式2.4.1，$\nabla\times(\nabla\phi)=0$，因此
    $$\nabla\times\boldsymbol{E} = -\nabla\times(\nabla\phi) - \nabla\times\frac{\partial \boldsymbol{A}}{\partial t} = -\frac{\partial}{\partial t}(\nabla\times\boldsymbol{A}) = -\frac{\partial \boldsymbol{B}}{\partial t}$$
    与法拉第电磁感应定律一致，证毕。

#### 4.3 四维时空与四维协变物理量
为构造洛伦兹不变的拉格朗日量，引入狭义相对论的四维时空形式，所有定义均与三维物理量严格对应。

**定义4.3.1 闵可夫斯基时空度规**：定义四维时空的度规张量为
$$\eta_{\mu\nu} = \eta^{\mu\nu} = \text{diag}(1, -1, -1, -1), \quad \mu,\nu=0,1,2,3$$
其中$\eta_{\mu\nu}$为协变度规，$\eta^{\mu\nu}$为逆变度规，满足$\eta^{\mu\lambda}\eta_{\lambda\nu} = \delta^\mu_\nu$，$\delta^\mu_\nu$为四维克罗内克δ函数（$\mu=\nu$时为1，否则为0）。

**定义4.3.2 四维时空坐标**：
- 逆变分量：$x^\mu = (x^0, x^1, x^2, x^3) = (ct, x, y, z)$
- 协变分量：$x_\mu = \eta_{\mu\nu}x^\nu = (ct, -x, -y, -z)$

**定义4.3.3 四维偏微分算子**：
- 协变分量：$\partial_\mu = \frac{\partial}{\partial x^\mu} = \left( \frac{\partial}{c\partial t}, \frac{\partial}{\partial x}, \frac{\partial}{\partial y}, \frac{\partial}{\partial z} \right) = \left( \frac{\partial}{c\partial t}, \nabla \right)$
- 逆变分量：$\partial^\mu = \eta^{\mu\nu}\partial_\nu = \left( \frac{\partial}{c\partial t}, -\nabla \right)$

**定义4.3.4 四维电磁势**：
- 逆变分量：$A^\mu = (A^0, A^1, A^2, A^3) = \left( \frac{\phi}{c}, A_x, A_y, A_z \right)$
- 协变分量：$A_\mu = \eta_{\mu\nu}A^\nu = \left( \frac{\phi}{c}, -A_x, -A_y, -A_z \right)$

**定义4.3.5 四维电流密度**：
- 逆变分量：$J^\mu = (J^0, J^1, J^2, J^3) = (c\rho, J_x, J_y, J_z)$
- 协变分量：$J_\mu = \eta_{\mu\nu}J^\nu = (c\rho, -J_x, -J_y, -J_z)$

**定理4.3.1 电荷守恒定律的四维形式**：$\partial_\mu J^\mu = 0$
**证明**：$\partial_\mu J^\mu = \partial_0 J^0 + \partial_1 J^1 + \partial_2 J^2 + \partial_3 J^3 = \frac{\partial}{c\partial t}(c\rho) + \nabla\cdot\boldsymbol{J} = \frac{\partial \rho}{\partial t} + \nabla\cdot\boldsymbol{J} = 0$，与三维连续性方程完全一致，证毕。

#### 4.4 电磁场张量（法拉第张量）
**定义4.4.1 电磁场张量**：由四维电磁势定义的二阶反对称张量，定义为
$$F_{\mu\nu} \equiv \partial_\mu A_\nu - \partial_\nu A_\mu \tag{4.4.1}$$
**性质4.4.1 反对称性**：$F_{\mu\nu} = -F_{\nu\mu}$，因此对角元$F_{\mu\mu}=0$（不对$\mu$求和）。

**定理4.4.1 电磁场张量与三维电场、磁场的对应关系**
电磁场张量的矩阵形式为：
$$F_{\mu\nu} = \begin{pmatrix}
0 & \frac{E_x}{c} & \frac{E_y}{c} & \frac{E_z}{c} \\
-\frac{E_x}{c} & 0 & -B_z & B_y \\
-\frac{E_y}{c} & B_z & 0 & -B_x \\
-\frac{E_z}{c} & -B_y & B_x & 0
\end{pmatrix} \tag{4.4.2}$$
**证明**：逐分量严格计算：
1.  时间-空间分量（$\mu=0, \nu=i, i=1,2,3$）：
    $$F_{0i} = \partial_0 A_i - \partial_i A_0$$
    代入$\partial_0=\frac{\partial}{c\partial t}$，$A_i = -A_{x_i}$，$A_0=\frac{\phi}{c}$，得：
    $$F_{0i} = \frac{\partial}{c\partial t}(-A_{x_i}) - \frac{\partial}{\partial x_i}\left( \frac{\phi}{c} \right) = -\frac{1}{c}\left( \frac{\partial A_{x_i}}{\partial t} + \frac{\partial \phi}{\partial x_i} \right)$$
    对比(4.2.1)式$E_i = -\left( \frac{\partial \phi}{\partial x_i} + \frac{\partial A_{x_i}}{\partial t} \right)$，因此$\frac{\partial A_{x_i}}{\partial t} + \frac{\partial \phi}{\partial x_i} = -E_i$，代入得：
    $$F_{0i} = -\frac{1}{c}(-E_i) = \frac{E_i}{c}$$
    由反对称性，$F_{i0} = -F_{0i} = -\frac{E_i}{c}$。

2.  空间-空间分量（$\mu=i, \nu=j, i,j=1,2,3, i\neq j$）：
    $$F_{ij} = \partial_i A_j - \partial_j A_i$$
    代入$A_j=-A_{x_j}, A_i=-A_{x_i}$，得：
    $$F_{ij} = \partial_i(-A_{x_j}) - \partial_j(-A_{x_i}) = -\left( \partial_i A_{x_j} - \partial_j A_{x_i} \right)$$
    由叉积的分量形式，$(\nabla\times\boldsymbol{A})_k = \sum_{i,j=1}^3 \varepsilon_{ijk} \partial_i A_{x_j}$，即$\partial_i A_{x_j} - \partial_j A_{x_i} = \varepsilon_{ijk} B_k$（$\boldsymbol{B}=\nabla\times\boldsymbol{A}$），因此：
    $$F_{ij} = -\varepsilon_{ijk} B_k$$
    例：$F_{12}=-\varepsilon_{123}B_3=-B_z$，$F_{13}=-\varepsilon_{132}B_2=B_y$，$F_{23}=-\varepsilon_{231}B_1=-B_x$，与矩阵(4.4.2)完全一致。

所有对角元$F_{\mu\mu}=\partial_\mu A_\mu - \partial_\mu A_\mu=0$，证毕。

**定义4.4.2 逆变电磁场张量**：$F^{\mu\nu} = \eta^{\mu\alpha}\eta^{\nu\beta}F_{\alpha\beta}$，其矩阵形式为：
$$F^{\mu\nu} = \begin{pmatrix}
0 & -\frac{E_x}{c} & -\frac{E_y}{c} & -\frac{E_z}{c} \\
\frac{E_x}{c} & 0 & -B_z & B_y \\
\frac{E_y}{c} & B_z & 0 & -B_x \\
\frac{E_z}{c} & -B_y & B_x & 0
\end{pmatrix} \tag{4.4.3}$$

**定理4.4.2 电磁场张量的洛伦兹标量**：
$$F_{\mu\nu}F^{\mu\nu} = 2\left( B^2 - \frac{E^2}{c^2} \right) \tag{4.4.4}$$
**证明**：
由爱因斯坦求和约定，$F_{\mu\nu}F^{\mu\nu} = \sum_{\mu=0}^3 \sum_{\nu=0}^3 F_{\mu\nu}F^{\mu\nu}$，利用反对称性，仅需计算$\mu<\nu$的项，再乘以2：
1.  时间-空间项（$\mu=0, \nu=i$）：$F_{0i}F^{0i} = (\frac{E_i}{c})(-\frac{E_i}{c}) = -\frac{E_i^2}{c^2}$，求和得$\sum_{i=1}^3 F_{0i}F^{0i} = -\frac{E^2}{c^2}$；
2.  空间-空间项（$\mu=i, \nu=j, i<j$）：$F_{12}F^{12}=(-B_z)(-B_z)=B_z^2$，$F_{13}F^{13}=(B_y)(B_y)=B_y^2$，$F_{23}F^{23}=(-B_x)(-B_x)=B_x^2$，求和得$B^2$。

因此总和为$2\times(-\frac{E^2}{c^2} + B^2) = 2(B^2 - \frac{E^2}{c^2})$，证毕。

---

### 第五章 电磁场的拉格朗日密度
#### 5.1 拉格朗日密度的构造原则
构造的电磁场拉格朗日密度必须满足三个核心原则：
1.  **洛伦兹不变性**：拉格朗日密度必须是洛伦兹标量，保证运动方程在所有惯性系中形式一致；
2.  **规范不变性**：在规范变换$A_\mu \to A_\mu + \partial_\mu \chi$（$\chi$为任意标量函数）下，拉格朗日密度的变分为四维散度，不影响运动方程；
3.  **线性场方程**：拉格朗日密度必须是场量的二阶项，保证导出的运动方程是线性的，符合电磁场的叠加原理。

#### 5.2 电磁场的总拉格朗日密度
总拉格朗日密度由自由电磁场项和场与电荷电流的相互作用项两部分组成：
$$\mathcal{L} = \mathcal{L}_0 + \mathcal{L}_{\text{int}} \tag{5.2.1}$$

##### 5.2.1 自由电磁场拉格朗日密度
满足上述原则的最低阶非平凡洛伦兹标量为$F_{\mu\nu}F^{\mu\nu}$，因此定义自由场拉格朗日密度为：
$$\mathcal{L}_0 = -\frac{1}{4\mu_0} F_{\mu\nu}F^{\mu\nu} \tag{5.2.2}$$
其中$\mu_0$为真空磁导率，系数$-\frac{1}{4\mu_0}$为归一化系数，保证导出的方程与实验规律一致。

将(4.4.4)式代入，利用真空中光速与介电常数、磁导率的关系$c = \frac{1}{\sqrt{\varepsilon_0\mu_0}}$（即$\frac{1}{c^2} = \varepsilon_0\mu_0$），可得三维形式：
$$\mathcal{L}_0 = -\frac{1}{4\mu_0} \cdot 2\left( B^2 - \frac{E^2}{c^2} \right) = \frac{1}{2}\varepsilon_0 E^2 - \frac{1}{2\mu_0} B^2 \tag{5.2.3}$$
该式与电磁场的能量密度经典形式完全一致，验证了构造的合理性。

##### 5.2.2 相互作用拉格朗日密度
场与电荷电流的相互作用项必须是洛伦兹标量，且与点电荷的洛伦兹力拉格朗日量一致，因此定义：
$$\mathcal{L}_{\text{int}} = - J^\mu A_\mu \tag{5.2.4}$$
该式在规范变换下，$\mathcal{L}_{\text{int}} \to -J^\mu (A_\mu + \partial_\mu \chi) = -J^\mu A_\mu - J^\mu \partial_\mu \chi$，第二项为四维散度$-\partial_\mu(J^\mu \chi) + \chi \partial_\mu J^\mu$，由电荷守恒$\partial_\mu J^\mu=0$，因此仅余四维散度项，对作用量的贡献为边界项，不影响运动方程，满足规范不变性。

---

## 第三篇 麦克斯韦方程组微分形式的完整推导
### 第六章 非齐次麦克斯韦方程组的推导
非齐次麦克斯韦方程组由场的欧拉-拉格朗日方程导出，每一步完全展开，无任何跳跃。

电磁场的独立场量为四维电磁势$A_\nu$（$\nu=0,1,2,3$），共4个独立分量，因此对每个分量$A_\nu$，都对应一个欧拉-拉格朗日方程：
$$\partial_\mu \left( \frac{\partial \mathcal{L}}{\partial (\partial_\mu A_\nu)} \right) - \frac{\partial \mathcal{L}}{\partial A_\nu} = 0, \quad \nu=0,1,2,3 \tag{6.0.1}$$

我们分两步计算方程中的两个偏导数项。

#### 步骤1：计算$\frac{\partial \mathcal{L}}{\partial A_\nu}$
总拉格朗日密度$\mathcal{L} = \mathcal{L}_0 + \mathcal{L}_{\text{int}}$中，仅相互作用项$\mathcal{L}_{\text{int}}=-J^\mu A_\mu$显含$A_\nu$，自由场项$\mathcal{L}_0$仅由$F_{\mu\nu}$（即$\partial_\mu A_\nu$）构成，不含$A_\nu$本身，因此：
$$\frac{\partial \mathcal{L}}{\partial A_\nu} = \frac{\partial}{\partial A_\nu} (-J^\mu A_\mu)$$
由爱因斯坦求和约定，$J^\mu A_\mu = \sum_{\mu=0}^3 J^\mu A_\mu$，对$A_\nu$求偏导，由克罗内克δ函数的筛选性：
$$\frac{\partial A_\mu}{\partial A_\nu} = \delta^\nu_\mu$$
因此：
$$\frac{\partial \mathcal{L}}{\partial A_\nu} = - J^\mu \delta^\nu_\mu = - J^\nu \tag{6.1.1}$$

#### 步骤2：计算$\frac{\partial \mathcal{L}}{\partial (\partial_\mu A_\nu)}$
仅自由场项$\mathcal{L}_0$含$\partial_\mu A_\nu$，因此：
$$\frac{\partial \mathcal{L}}{\partial (\partial_\mu A_\nu)} = \frac{\partial \mathcal{L}_0}{\partial (\partial_\mu A_\nu)} = -\frac{1}{4\mu_0} \cdot \frac{\partial}{\partial (\partial_\mu A_\nu)} \left( F_{\alpha\beta}F^{\alpha\beta} \right) \tag{6.2.1}$$

分步计算偏导数：
##### 子步骤2.1 计算$\frac{\partial F_{\alpha\beta}}{\partial (\partial_\mu A_\nu)}$
由电磁场张量的定义$F_{\alpha\beta} = \partial_\alpha A_\beta - \partial_\beta A_\alpha$，对$\partial_\mu A_\nu$求偏导：
$$\frac{\partial F_{\alpha\beta}}{\partial (\partial_\mu A_\nu)} = \frac{\partial}{\partial (\partial_\mu A_\nu)}(\partial_\alpha A_\beta) - \frac{\partial}{\partial (\partial_\mu A_\nu)}(\partial_\beta A_\alpha)$$
偏导数的筛选性：$\frac{\partial (\partial_\alpha A_\beta)}{\partial (\partial_\mu A_\nu)} = \delta^\mu_\alpha \delta^\nu_\beta$（仅当$\alpha=\mu$且$\beta=\nu$时，偏导数为1，否则为0），因此：
$$\frac{\partial F_{\alpha\beta}}{\partial (\partial_\mu A_\nu)} = \delta^\mu_\alpha \delta^\nu_\beta - \delta^\mu_\beta \delta^\nu_\alpha \tag{6.2.2}$$

##### 子步骤2.2 计算$\frac{\partial (F_{\alpha\beta}F^{\alpha\beta})}{\partial (\partial_\mu A_\nu)}$
由乘积求导法则（高数内容）：
$$\frac{\partial (F_{\alpha\beta}F^{\alpha\beta})}{\partial (\partial_\mu A_\nu)} = F_{\alpha\beta} \frac{\partial F^{\alpha\beta}}{\partial (\partial_\mu A_\nu)} + F^{\alpha\beta} \frac{\partial F_{\alpha\beta}}{\partial (\partial_\mu A_\nu)}$$
由于$F^{\alpha\beta}$与$F_{\alpha\beta}$仅差度规张量的常数收缩，二者对$\partial_\mu A_\nu$的偏导数具有相同的线性关系，因此两项相等，即：
$$\frac{\partial (F_{\alpha\beta}F^{\alpha\beta})}{\partial (\partial_\mu A_\nu)} = 2 F^{\alpha\beta} \frac{\partial F_{\alpha\beta}}{\partial (\partial_\mu A_\nu)} \tag{6.2.3}$$

将(6.2.2)式代入(6.2.3)式：
$$\frac{\partial (F_{\alpha\beta}F^{\alpha\beta})}{\partial (\partial_\mu A_\nu)} = 2 F^{\alpha\beta} \left( \delta^\mu_\alpha \delta^\nu_\beta - \delta^\mu_\beta \delta^\nu_\alpha \right)$$
由克罗内克δ函数的筛选性，第一项中$\alpha=\mu, \beta=\nu$，得$F^{\mu\nu}$；第二项中$\beta=\mu, \alpha=\nu$，得$F^{\nu\mu}$，因此：
$$\frac{\partial (F_{\alpha\beta}F^{\alpha\beta})}{\partial (\partial_\mu A_\nu)} = 2 F^{\mu\nu} - 2 F^{\nu\mu}$$
由电磁场张量的反对称性$F^{\nu\mu} = -F^{\mu\nu}$，代入得：
$$\frac{\partial (F_{\alpha\beta}F^{\alpha\beta})}{\partial (\partial_\mu A_\nu)} = 2 F^{\mu\nu} - 2 (-F^{\mu\nu}) = 4 F^{\mu\nu} \tag{6.2.4}$$

##### 子步骤2.3 代入得到最终结果
将(6.2.4)式代入(6.2.1)式：
$$\frac{\partial \mathcal{L}}{\partial (\partial_\mu A_\nu)} = -\frac{1}{4\mu_0} \cdot 4 F^{\mu\nu} = -\frac{F^{\mu\nu}}{\mu_0} \tag{6.2.5}$$

#### 步骤3：代入欧拉-拉格朗日方程
将(6.1.1)式和(6.2.5)式代入欧拉-拉格朗日方程(6.0.1)：
$$\partial_\mu \left( -\frac{F^{\mu\nu}}{\mu_0} \right) - (-J^\nu) = 0$$
整理得：
$$-\frac{1}{\mu_0} \partial_\mu F^{\mu\nu} + J^\nu = 0$$
两边乘以$-\mu_0$，得到麦克斯韦方程组的**四维协变非齐次形式**：
$$\partial_\mu F^{\mu\nu} = \mu_0 J^\nu \tag{6.3.1}$$

#### 步骤4：将四维方程分解为三维微分形式
四维方程(6.3.1)对每个自由指标$\nu=0,1,2,3$都成立，我们分别取$\nu=0$（时间分量）和$\nu=i$（空间分量），分解为两个三维方程。

##### 6.4.1 高斯定理的微分形式
取$\nu=0$，代入(6.3.1)式，由爱因斯坦求和约定，对$\mu=0,1,2,3$求和：
$$\partial_0 F^{00} + \partial_1 F^{10} + \partial_2 F^{20} + \partial_3 F^{30} = \mu_0 J^0$$
由电磁场张量的反对称性，$F^{00}=0$，因此第一项为0；
代入$F^{i0} = \frac{E_i}{c}$（由修正后的逆变张量矩阵），$J^0 = c\rho$，得：
$$\partial_1 \left( \frac{E_x}{c} \right) + \partial_2 \left( \frac{E_y}{c} \right) + \partial_3 \left( \frac{E_z}{c} \right) = \mu_0 \cdot c\rho$$
左边提取公因子$\frac{1}{c}$，得：
$$\frac{1}{c} \left( \frac{\partial E_x}{\partial x} + \frac{\partial E_y}{\partial y} + \frac{\partial E_z}{\partial z} \right) = \mu_0 c \rho$$
左边括号内为$\nabla\cdot\boldsymbol{E}$，因此：
$$\frac{1}{c} \nabla\cdot\boldsymbol{E} = \mu_0 c \rho$$
两边乘以$c$，得：
$$\nabla\cdot\boldsymbol{E} = \mu_0 c^2 \rho$$
代入$c^2 = \frac{1}{\varepsilon_0 \mu_0}$，得$\mu_0 c^2 = \frac{1}{\varepsilon_0}$，因此得到**高斯定理的微分形式**：
$$\nabla\cdot\boldsymbol{E} = \frac{\rho}{\varepsilon_0} \tag{6.4.1}$$

##### 6.4.2 安培-麦克斯韦环路定理的微分形式
取$\nu=i$（$i=1,2,3$，空间分量），代入四维方程(6.3.1)：
$$\partial_\mu F^{\mu i} = \mu_0 J^i$$
对$\mu=0,1,2,3$求和，拆分为时间分量$\mu=0$和空间分量$\mu=j$（$j=1,2,3$）：
$$\partial_0 F^{0i} + \sum_{j=1}^3 \partial_j F^{ji} = \mu_0 J^i \tag{6.4.2}$$

分别计算两项：
1.  第一项：$\partial_0 F^{0i}$
    代入$\partial_0 = \frac{\partial}{c\partial t}$，$F^{0i} = -\frac{E_i}{c}$，得：
    $$\partial_0 F^{0i} = \frac{\partial}{c\partial t} \left( -\frac{E_i}{c} \right) = -\frac{1}{c^2} \frac{\partial E_i}{\partial t} \tag{6.4.3}$$

2.  第二项：$\sum_{j=1}^3 \partial_j F^{ji}$
    由逆变张量矩阵，空间-空间分量$F^{ji} = -\varepsilon_{jik} B_k$，代入得：
    $$\sum_{j=1}^3 \partial_j F^{ji} = -\sum_{j=1}^3 \partial_j (\varepsilon_{jik} B_k) = -\sum_{j,k=1}^3 \varepsilon_{jik} \partial_j B_k$$
    由列维-奇维塔符号的反对称性$\varepsilon_{jik} = -\varepsilon_{ijk}$，代入得：
    $$\sum_{j=1}^3 \partial_j F^{ji} = \sum_{j,k=1}^3 \varepsilon_{ijk} \partial_j B_k = (\nabla\times\boldsymbol{B})_i \tag{6.4.4}$$

将(6.4.3)和(6.4.4)代入(6.4.2)式，得：
$$-\frac{1}{c^2} \frac{\partial E_i}{\partial t} + (\nabla\times\boldsymbol{B})_i = \mu_0 J_i$$
移项得：
$$(\nabla\times\boldsymbol{B})_i = \mu_0 J_i + \frac{1}{c^2} \frac{\partial E_i}{\partial t}$$
代入$\frac{1}{c^2} = \varepsilon_0 \mu_0$，写成矢量形式，即**安培-麦克斯韦环路定理的微分形式**：
$$\nabla\times\boldsymbol{B} = \mu_0 \boldsymbol{J} + \mu_0 \varepsilon_0 \frac{\partial \boldsymbol{E}}{\partial t} \tag{6.4.5}$$

---

### 第七章 齐次麦克斯韦方程组的推导
齐次麦克斯韦方程组直接由电磁场张量的定义导出，本质是电磁场的Bianchi恒等式。

#### 7.1 Bianchi恒等式的证明
**定理7.1.1 电磁场张量的Bianchi恒等式**
对任意指标$\lambda,\mu,\nu$，有
$$\partial_\lambda F_{\mu\nu} + \partial_\mu F_{\nu\lambda} + \partial_\nu F_{\lambda\mu} = 0 \tag{7.1.1}$$
**证明**：
将电磁场张量的定义$F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$代入左边，完全展开：
$$\partial_\lambda (\partial_\mu A_\nu - \partial_\nu A_\mu) + \partial_\mu (\partial_\nu A_\lambda - \partial_\lambda A_\nu) + \partial_\nu (\partial_\lambda A_\mu - \partial_\mu A_\lambda)$$
展开后共6项：
$$\partial_\lambda\partial_\mu A_\nu - \partial_\lambda\partial_\nu A_\mu + \partial_\mu\partial_\nu A_\lambda - \partial_\mu\partial_\lambda A_\nu + \partial_\nu\partial_\lambda A_\mu - \partial_\nu\partial_\mu A_\lambda \tag{7.1.2}$$

由高等数学中的**混合偏导数相等定理**：若函数$A_\nu$的二阶偏导数连续，则$\partial_\lambda\partial_\mu A_\nu = \partial_\mu\partial_\lambda A_\nu$，即偏导数的顺序可交换。

将(7.1.2)式的项两两配对：
1.  第一项$\partial_\lambda\partial_\mu A_\nu$与第四项$-\partial_\mu\partial_\lambda A_\nu$：由混合偏导数相等，二者相等，相加为0；
2.  第二项$-\partial_\lambda\partial_\nu A_\mu$与第五项$\partial_\nu\partial_\lambda A_\mu$：同理，相加为0；
3.  第三项$\partial_\mu\partial_\nu A_\lambda$与第六项$-\partial_\nu\partial_\mu A_\lambda$：同理，相加为0。

因此(7.1.2)式的总和为0，Bianchi恒等式得证。

#### 7.2 将Bianchi恒等式分解为三维微分形式
我们将Bianchi恒等式按指标的取值分为两种情况，分别对应两个齐次麦克斯韦方程。

##### 7.2.1 法拉第电磁感应定律的微分形式
取Bianchi恒等式的指标为：一个时间指标，两个空间指标，即$\lambda=0, \mu=i, \nu=j$（$i,j=1,2,3$），代入(7.1.1)式：
$$\partial_0 F_{ij} + \partial_i F_{j0} + \partial_j F_{0i} = 0 \tag{7.2.1}$$

代入各分量的表达式：
1.  $F_{ij} = -\varepsilon_{ijk} B_k$；
2.  $F_{j0} = -F_{0j} = -\frac{E_j}{c}$；
3.  $F_{0i} = \frac{E_i}{c}$。

将以上代入(7.2.1)式，两边乘以$c$消去分母，得：
$$\partial_0 (-\varepsilon_{ijk} B_k) \cdot c + \partial_i (-E_j) + \partial_j E_i = 0$$
代入$\partial_0 = \frac{\partial}{c\partial t}$，$\partial_0 \cdot c = \frac{\partial}{\partial t}$，化简得：
$$-\varepsilon_{ijk} \frac{\partial B_k}{\partial t} + \partial_j E_i - \partial_i E_j = 0$$

取具体指标验证：令$i=1,j=2$，则$\varepsilon_{123}=1$，代入得：
$$-\frac{\partial B_z}{\partial t} + \frac{\partial E_x}{\partial y} - \frac{\partial E_y}{\partial x} = 0$$
整理得：
$$\frac{\partial E_y}{\partial x} - \frac{\partial E_x}{\partial y} = -\frac{\partial B_z}{\partial t}$$
而左边正是$(\nabla\times\boldsymbol{E})_z$，因此$(\nabla\times\boldsymbol{E})_z = -\frac{\partial B_z}{\partial t}$。

同理，取其余指标组合，可得到$x,y$分量的对应关系，写成矢量形式，即**法拉第电磁感应定律的微分形式**：
$$\nabla\times\boldsymbol{E} = -\frac{\partial \boldsymbol{B}}{\partial t} \tag{7.2.2}$$

##### 7.2.2 磁场高斯定理的微分形式
取Bianchi恒等式的指标为全空间指标，即$\lambda=i, \mu=j, \nu=k$（$i,j,k=1,2,3$），代入(7.1.1)式：
$$\partial_i F_{jk} + \partial_j F_{ki} + \partial_k F_{ij} = 0 \tag{7.2.3}$$

取具体指标$i=1,j=2,k=3$（其余置换结果一致），代入得：
$$\partial_1 F_{23} + \partial_2 F_{31} + \partial_3 F_{12} = 0$$
代入$F_{23}=-B_x$，$F_{31}=-B_y$，$F_{12}=-B_z$，得：
$$\partial_1 (-B_x) + \partial_2 (-B_y) + \partial_3 (-B_z) = 0$$
两边乘以-1，得：
$$\frac{\partial B_x}{\partial x} + \frac{\partial B_y}{\partial y} + \frac{\partial B_z}{\partial z} = 0$$
写成矢量形式，即**磁场高斯定理的微分形式**：
$$\nabla\cdot\boldsymbol{B} = 0 \tag{7.2.4}$$

---

#### 麦克斯韦方程组微分形式汇总
至此，我们完整推导出真空中麦克斯韦方程组的全部四个微分形式：
1.  **电场高斯定理**：$\displaystyle \nabla\cdot\boldsymbol{E} = \frac{\rho}{\varepsilon_0}$
2.  **磁场高斯定理**：$\displaystyle \nabla\cdot\boldsymbol{B} = 0$
3.  **法拉第电磁感应定律**：$\displaystyle \nabla\times\boldsymbol{E} = -\frac{\partial \boldsymbol{B}}{\partial t}$
4.  **安培-麦克斯韦环路定理**：$\displaystyle \nabla\times\boldsymbol{B} = \mu_0\boldsymbol{J} + \mu_0\varepsilon_0\frac{\partial \boldsymbol{E}}{\partial t}$

---

## 第四篇 麦克斯韦方程组积分形式的完整推导
麦克斯韦方程组的积分形式由微分形式结合矢量分析的高斯散度定理和斯托克斯旋度定理导出，每一步均与高数内容严格对应。

### 第八章 积分形式的逐式推导
#### 8.1 电场高斯定理的积分形式
**微分形式**：$\nabla\cdot\boldsymbol{E} = \frac{\rho}{\varepsilon_0}$

在三维空间中取任意有界闭区域$V$，其边界为分片光滑的闭合曲面$S$，对微分形式两边在$V$上做体积分：
$$\int_V \nabla\cdot\boldsymbol{E} \, dV = \int_V \frac{\rho}{\varepsilon_0} dV \tag{8.1.1}$$

1.  左边：由**高斯散度定理**，$\int_V \nabla\cdot\boldsymbol{E} dV = \oint_S \boldsymbol{E}\cdot d\boldsymbol{S}$，即电场强度$\boldsymbol{E}$在闭合曲面$S$上的电通量；
2.  右边：提取常数$1/\varepsilon_0$，$\int_V \rho dV = Q$，即闭合曲面$S$所包围的总电荷量。

代入(8.1.1)式，得**电场高斯定理的积分形式**：
$$\oint_S \boldsymbol{E}\cdot d\boldsymbol{S} = \frac{Q}{\varepsilon_0} \tag{8.1.2}$$

#### 8.2 磁场高斯定理的积分形式
**微分形式**：$\nabla\cdot\boldsymbol{B} = 0$

对任意闭合曲面$S$包围的体积$V$，两边做体积分：
$$\int_V \nabla\cdot\boldsymbol{B} \, dV = 0$$

由高斯散度定理，左边为磁感应强度$\boldsymbol{B}$在闭合曲面$S$上的磁通量，因此得**磁场高斯定理的积分形式**：
$$\oint_S \boldsymbol{B}\cdot d\boldsymbol{S} = 0 \tag{8.2.1}$$

#### 8.3 法拉第电磁感应定律的积分形式
**微分形式**：$\nabla\times\boldsymbol{E} = -\frac{\partial \boldsymbol{B}}{\partial t}$

取任意分片光滑的有向曲面$\Sigma$，其边界为分段光滑的闭合有向曲线$L$，$L$的正向与$\Sigma$的法向满足右手螺旋定则，对微分形式两边在$\Sigma$上做曲面积分：
$$\int_\Sigma (\nabla\times\boldsymbol{E})\cdot d\boldsymbol{S} = \int_\Sigma \left( -\frac{\partial \boldsymbol{B}}{\partial t} \right) \cdot d\boldsymbol{S} \tag{8.3.1}$$

1.  左边：由**斯托克斯旋度定理**，$\int_\Sigma (\nabla\times\boldsymbol{E})\cdot d\boldsymbol{S} = \oint_L \boldsymbol{E}\cdot d\boldsymbol{l}$，即电场强度$\boldsymbol{E}$沿闭合曲线$L$的环量；
2.  右边：积分区域$\Sigma$不随时间变化，偏导数与积分可交换顺序（高数含参变量积分的可微性定理），即：
    $$\int_\Sigma \frac{\partial \boldsymbol{B}}{\partial t} \cdot d\boldsymbol{S} = \frac{d}{dt} \int_\Sigma \boldsymbol{B}\cdot d\boldsymbol{S}$$
    其中$\Phi_m = \int_\Sigma \boldsymbol{B}\cdot d\boldsymbol{S}$为穿过曲面$\Sigma$的磁通量。

代入(8.3.1)式，得**法拉第电磁感应定律的积分形式**：
$$\oint_L \boldsymbol{E}\cdot d\boldsymbol{l} = -\frac{d\Phi_m}{dt} \tag{8.3.2}$$

#### 8.4 安培-麦克斯韦环路定理的积分形式
**微分形式**：$\nabla\times\boldsymbol{B} = \mu_0 \boldsymbol{J} + \mu_0 \varepsilon_0 \frac{\partial \boldsymbol{E}}{\partial t}$

对以闭合曲线$L$为边界的任意曲面$\Sigma$，两边做曲面积分：
$$\int_\Sigma (\nabla\times\boldsymbol{B})\cdot d\boldsymbol{S} = \int_\Sigma \left( \mu_0 \boldsymbol{J} + \mu_0 \varepsilon_0 \frac{\partial \boldsymbol{E}}{\partial t} \right) \cdot d\boldsymbol{S} \tag{8.4.1}$$

1.  左边：由斯托克斯旋度定理，$\int_\Sigma (\nabla\times\boldsymbol{B})\cdot d\boldsymbol{S} = \oint_L \boldsymbol{B}\cdot d\boldsymbol{l}$，即磁感应强度$\boldsymbol{B}$沿闭合曲线$L$的环量；
2.  右边第一项：$\int_\Sigma \boldsymbol{J}\cdot d\boldsymbol{S} = I$，即穿过曲面$\Sigma$的传导电流；
3.  右边第二项：提取常数$\mu_0 \varepsilon_0$，交换偏导数与积分顺序，得$\mu_0 \varepsilon_0 \frac{d}{dt} \int_\Sigma \boldsymbol{E}\cdot d\boldsymbol{S} = \mu_0 \varepsilon_0 \frac{d\Phi_e}{dt}$，其中$\Phi_e = \int_\Sigma \boldsymbol{E}\cdot d\boldsymbol{S}$为穿过曲面$\Sigma$的电通量，$\varepsilon_0 \frac{d\Phi_e}{dt}$为位移电流。

代入(8.4.1)式，得**安培-麦克斯韦环路定理的积分形式**：
$$\oint_L \boldsymbol{B}\cdot d\boldsymbol{l} = \mu_0 \left( I + \varepsilon_0 \frac{d\Phi_e}{dt} \right) \tag{8.4.2}$$

---

#### 麦克斯韦方程组积分形式汇总
1.  **电场高斯定理**：$\displaystyle \oint_S \boldsymbol{E}\cdot d\boldsymbol{S} = \frac{Q}{\varepsilon_0}$
2.  **磁场高斯定理**：$\displaystyle \oint_S \boldsymbol{B}\cdot d\boldsymbol{S} = 0$
3.  **法拉第电磁感应定律**：$\displaystyle \oint_L \boldsymbol{E}\cdot d\boldsymbol{l} = -\frac{d\Phi_m}{dt}$
4.  **安培-麦克斯韦环路定理**：$\displaystyle \oint_L \boldsymbol{B}\cdot d\boldsymbol{l} = \mu_0 \left( I + \varepsilon_0 \frac{d\Phi_e}{dt} \right)$
