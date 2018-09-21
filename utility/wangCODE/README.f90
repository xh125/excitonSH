!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% this program is used to simulate charge transport with the self-consistent surface hopping method %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% one-dimensional molecular stack with local electron-phonon couplings and system-bath interactions %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% by the wang group at department of chemistry, zhejiang university; 2017/03/09; ljwang@zju.edu.cn  %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

the input file name is 'SHIN'

!变量的定义：
integer
dq=kind(1.0D0) !双精度Kind的设置值
dqc=kind((1.0D0,1.0D0)) !双精度复数的kind的设置

na1site,	!二维(1 and 3)情况下在a1方向上的超胞中原胞个数
ia1site
na2site   !二维(1 and 3)情况下在a2方向上的超胞中原胞个数,1D情况下，na3site=1，na2site=1
ia2site		
na3site		!二维(1 and 3)情况下在a3方向上的超胞中原胞个数,2D情况下，na3site=1
ia3site

nwann,		!	一个原胞内的wannier基矢数目。
!nsite		!整个超胞中的基矢数目=nxsite*nysite*nzsite*nwann
nbasis		!整个超胞中的基矢数目=na1site*na2site*na3site*nwann
ibasis    !

natoms,		!一个原胞内的原子数
nfreem,		!一个原胞内原子核自由度=3*natoms
ifreem,

!nsite，   !一维原子链的原子数,!二维情况下的H的维度。(2nldim+1)*(2nhdim+1)*nwann*nwann
naver，		!一共计算naver条轨迹
nsnap，		!每一条轨迹的MD的帧数
nstep,    !在每一步的核的MD过程中进行nstep步电子的量子动力学计算
iaver,    !第iaver条轨迹
isnap,    !核MD的第isnap步
istep,		!核一次MD步长过程中第istep步电子步数

isurface, !  初始电子所在的势能面（绝热势能面），电子（空穴）随时间的占据情况。

!icenter,	!icenter=anint(nsite/2.0d0)  初始电子所在的site，用于描述初始电子态在透热表象下的状态
!icenter_a1,
!icenter_a2,
!icenter_a3,
icenter(4) !表示在绝热表象下的波函数，a1,a2,a3,iwann
icenter_index,
init_wann

!isite,		!第isite个site位置
!!use in read wannier90_hr.dat
nrpts,		!wannier计算临近单胞间的转移积分，临近的单胞个数（包含自己）
irpts,
innwan,		!转移积分矩阵右失
imnwan，	!转移积分矩阵左矢
n					!
m					!
ir1
ir2
ir3
ReH
ImH

real(kind=8) !原子坐标到各单位坐标的变换系数
au2cm,
au2ev,
au2j,
au2fs,
au2ps,
au2amu,
au2ang,
kb,				!玻尔兹曼常数 kb=1.3806504d-23
sqrt3,
sqrt5,
sqrt7

real(kind=8) 
!mass,			!每个site处核沿简正坐标振动的有效质量
!mass(:),		!mass(nfreem) 每个ifreem自由度沿简正坐标振动的有效质量
temp,			!体系温度T，用于随机产生初始速度和位置的高斯分布（玻尔兹曼分布）的展宽
gamma,		!体系受力情况的阻力系数-> dv(isite)=(-k*xx(isite)-alpha*pp(isite,isurface)**2)/mass-gamma*vv(isite)，&
          !单位为1/t
dt,				!量子动力学的时间步长,没两帧之间有nstep步电子动力学，步长为dt
!alpha,		!local电声耦合强度
!k,				!体系沿简正坐标振动的回复力系数，F=-kdQ
!k(:),			!k(nfreem)	体系ifreem自由度沿简正坐标振动的回复力系数，F=-kdQ
!tau				!最邻近基矢之间的转移积分，体现在H对角线两侧

Latt(3,3) !体系晶格矢量
Rwann(:,:) !Rwann(3,nwann) home单胞中wannier波函数的中心。
real(kind=8),allocatable :: 
h0(:,:),	!h0(nsite,nsite) 初始的H
hh(:,:)		!(hh(1:nsite,1:nsite)) 考虑电声耦合之后的H
hep(:,:,:)	!hep(nsite,nsite,nfreem) 电声耦合强度矩阵
real(kind=8),allocatable ::
 !x(:),		!x(nsite) 体系的位置
 !v(:),		!v(nsite) 体系的速度
 q(:)			!q(nfreem) 体系简正坐标位置
 v(:)			!v(nfreem) 体系简正坐标速度
 e(:),		!e(nsite) 体系的能量本征值
 p(:,:),	!p(nsite,nsite) 体系的能量本征态
 !d(:,:,:), !d(nsite,nsite,nsite) 体系不同本征态之间在不同速度方向上的非绝热耦合项。
 d(:,:,:) !dd(nsite,nsite,nfreem) 体系不同本征态沿某一自由度方向上的非绝热耦合项
 g(:)			!g(nsite)  gg(isite)=2.0d0*tt*real(conjg(ww(isurface))*ww(isite))*sumvd/real(conjg(ww(isurface))*ww(isurface))
					!跃迁几率
real(kind=8),allocatable ::
 !x0(:),   !x0(nsite)
 !v0(:),		!v0(nsite)
 q0(:)			!q0(nfreem)
 v0(:)			!v0(nfreem)
 e0(:),		!e0(nsite)
 p0(:,:),	!p0(nsite,nsite)
 d0(:,:,:), !d0(nsite,nsite,nsite) 非绝热耦合项   ???????
						!dd(isite,jsite,ksite)=alpha*pp(ksite,isite)*pp(ksite,jsite)/(ee(jsite)-ee(isite))
 g1(:)		!g1(nsite) 在绝热表象下的跃迁几率
real(kind=8),allocatable :: 
pes(:,:,:),	! pes(0:nsite,1:nsnap,1:naver） potential energy surface各个iaver轨线在MD过程中的各个势能面，即各个isnap中的能量本征值随MD的变换
						! pes(0,isnap,iaver)=e(isurface) 当前电子所在的势能面
						! pes(isite,isnap,iaver)=e(isite) 所有势能面
inf(:,:,:),	!	inf(3,nsnap,naver)
						!inf(1,isnap,iaver)=sumg0
						!inf(2,isnap,iaver)=sumg1
						!inf(3,isnap,iaver)=minde
csit(:,:),	!	csit(nsite,nsnap) 电子在nsnap步时，在nsite局域轨道(非绝热表象)下的平均占据数
wsit(:,:),	!	wsite(nsite,nsnap) 电子在nsnap步时，在nsite能量本征态(绝热表象)下的平均占据数
psit(:,:),	! psite(nsite,nsnap)
xsit(:,:),	!	xsite(nsite,nsnap)
ksit(:,:),	!	ksit(nsite,nsnap)
msd(:),			! msd(nsnap) 平均的扩散 r2
ipr(:),			!	ipr(nsnap)
msds(:,:)		!	msds(nsnap,naver) 每一条轨迹的<|r2|>
complex(kind=16) eye  !eye=(0.0d0,1.0d0) 复数i单位
complex(kind=16),allocatable ::
c(:),				!c(nsite) cc=0.0d0; cc(icenter)=1.0d0 在透热表象下的电子态占据表示
w(:),				!w(nsite)	在绝热表象（能量本征态）下的电子态	
w0(:)				!w0(nsite)	初始电子态在绝热表象下的展开系数

real(kind=8) 
t0,			!CPU的初始时间（以秒计算）
t1      !CPU的结束时间（以秒计算）用于计算计算所花时间
real(kind=8) 
sumg0,			!在非绝热表象下总的跃迁几率
sumg1,			!在绝热表象下总的跃迁几率
minde,			!临近的两个势能面之间的能量差
flagd				!用于存放随机数或者其他变量
flagr				!用于存放（0-1）之间均匀分布的随机数

real(kind=8) 
pp(1:nsite,1:nsite),
tt,  !->dt
tt2,	!->dt/2
tt6		!->dt/6
real(kind=8) 
xx(nsite),
xx0(nsite),
dx1(nsite),
dx2(nsite),
dx3(nsite),
dx4(nsite)
real(kind=8)
v(nsite),
vv0(nsite),
dv1(nsite),
dv2(nsite),
dv3(nsite),
dv4(nsite)

complex*16 
cc(1:nsite),
cc0(1:nsite),
dc1(1:nsite),
dc2(1:nsite),
dc3(1:nsite),
dc4(1:nsite)

real(kind=8)
sumvd					!

real(kind=8)
sigmar,				! sigmar=dsqrt(2.0d0*kb*temp*gamma*tt/mass)
kk,						! 
							!   kk=k
							!		do jsite=1,nsite
							!			if(jsite.ne.isurface) kk=kk+2.0d0*alpha*dd(jsite,isurface,isite)*pp(isite,isurface)*pp(isite,jsite)
							!		enddo
r1,						!
r2,						!
r3,						!
r4,						!
z1,						!
z2,						!
z3,						!
z4,						!
gaussian_random_number_fast		!
