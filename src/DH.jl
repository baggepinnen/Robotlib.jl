mutable struct DH{T}
    dhpar::Matrix{T}
    offset::Vector{T}
    GR::Matrix{T}
end

DH(dhpar::AbstractMatrix{T},offset) where T = DH(dhpar,offset,Matrix{T}(I, length(offset),length(offset)))

function DH2400()
    DH_alpha0 = -π/2
    DH_alpha1 =  0
    DH_alpha2 = -π/2
    DH_alpha3 = -π/2
    DH_alpha4 = +π/2
    DH_alpha5 =  0

    DH_a0= 100
    DH_a1= 705
    DH_a2= 135

    DH_d0 =615
    DH_d1 =0
    DH_d2 =0
    DH_d3 =755
    DH_d4 =0
    DH_d5 =85

    offset = [0, -90, 0, 180, 0, 0].*π/180

    dhpar   = [ DH_alpha0      DH_a0        0         DH_d0        0;
    DH_alpha1      DH_a1        0         DH_d1        0;
    DH_alpha2      DH_a2        0         DH_d2        0;
    DH_alpha3       0           0         DH_d3        0;
    DH_alpha4       0           0         DH_d4        0;
    DH_alpha5       0           0         DH_d5        0 ]

    dhpar[:,[2,4]] /= 1000
    return DH(dhpar,offset)
end

function DHYuMi(left=true)
    # 100 100 -100 100 -101 100 100
    # baseAnglesRight = [0.63 ; 0.95 ; 0.18 ];
    # rotRight = rotx(0.63)*roty(0.95)*rotz(0.18);
    # tRight = [39.3 -64.5 403.5]';
    # baseRight = [rotRight tRight; 0 0 0 1]
    #
    # baseAnglesLeft = [-0.629355 ; 0.950621 ; -0.184101 ];
    # rotLeft = rotx(-0.629355)*roty(0.950621)*rotz(-0.184101);
    # tLeft = [39.1 66.5 407.3]';
    # baseLeft = [rotLeft tLeft; 0 0 0 1]
    baseAnglesRight = [0.63 , 0.95 , 0.18 ]
    baseAnglesLeft = [-0.63 , 0.95 , -0.18]

    # gear ratio must be applied AFTER abb2logical
    gearRatio = zeros(7,7)
    gearRatio[1,1] = 1/100
    gearRatio[2,2] = 1/100

    gearRatio[7,7] = 1/100
    gearRatio[3,3] = 1/100
    gearRatio[4,4] = -1/100
    gearRatio[5,5] = 1/100
    gearRatio[6,6] = -1/101

#     tmp(1:2) = tmp(1:2)/100;
# tmp(3) = -1*tmp(3)/100;
# tmp(4) = tmp(4)/100;
# tmp(5) = -1*tmp(5)/101;
# tmp(6:7) = tmp(6:7)/100;

    DH_alpha0 = -pi/2;
    DH_alpha1 = pi/2;
    DH_alpha2 = -pi/2;
    DH_alpha3 = pi/2;
    DH_alpha4 = -pi/2;
    DH_alpha5 =  pi/2;
    DH_alpha6 = 0;

    upArmLen =      0.2465;
    lowArmLen =     0.265;
    handLen =       0.032;
    shoulderLen =   0.110;

    shoulderOffs = -0.030;
    elbowOffs =     0.0405;
    wristOffs =     0.0135;
    handOffs =     -0.027;
    n_joints = 7;
    offset = [0,0,0,pi/2,0,0,0]

    #          alpha                a           theta       d           Rev=0/Pris=1    theta-offset
    dhpar   = [ DH_alpha0      shoulderOffs     0       shoulderLen         0
                DH_alpha1     -shoulderOffs     0           0               0
                DH_alpha2      elbowOffs        0       upArmLen            0
                DH_alpha3     -elbowOffs        0           0               0
                DH_alpha4      wristOffs        0       lowArmLen           0
                DH_alpha5      handOffs         0           0               0
                DH_alpha6      0                0       handLen             0  ];

    return DH(dhpar,offset,gearRatio)
end


function DH7600()
    DH_alpha= [-pi/2 0 -pi/2 pi/2 -pi/2 0];
    DH_a0= 410;
    DH_a1= 1075;
    DH_a2= 165;

    DH_d0 =780;
    DH_d1 =0;
    DH_d2 =0;
    DH_d3 = 1056;
    DH_d4 =0;
    DH_d5 =250;

    offset = [0, -90, 0, 0, 0, 180]*pi/180;

    dhpar   = [ DH_alpha[1]      DH_a0        0         DH_d0        0
    DH_alpha[2]      DH_a1        0         DH_d1        0
    DH_alpha[3]      DH_a2        0         DH_d2        0
    DH_alpha[4]       0           0         DH_d3        0
    DH_alpha[5]       0           0         DH_d4        0
    DH_alpha[6]       0           0         DH_d5        0 ];

    dhpar[:,[2,4]] /= 1000

    return DH(dhpar,offset)
end

function DHtest()
    DH_alpha = [1.0,2,3,4,5,6]
    DH_a0    = 0;
    DH_a1    = 0;
    DH_a2    = 0;

    DH_d0    = 780;
    DH_d1    = 10;
    DH_d2    = 20;
    DH_d3    = 1056;
    DH_d4    = 30;
    DH_d5    = 250;

    offset = zeros(6);

    dhpar   = [ DH_alpha[1]      DH_a0        0         DH_d0        0
    DH_alpha[2]      DH_a1        0         DH_d1        0
    DH_alpha[3]      DH_a2        0         DH_d2        0
    DH_alpha[4]       0           0         DH_d3        0
    DH_alpha[5]       0           0         DH_d4        0
    DH_alpha[6]       0           0         DH_d5        0 ];

    return DH(dhpar/1000,offset)
end


function abb2logical!(q::AbstractMatrix)
    q .= q[:,SA[1, 2, 7, 3, 4, 5, 6]]
end

function logical2abb!(q::AbstractMatrix)
    q .= q[:,SA[1, 2, 4, 5, 6, 7, 3]]
end

function abb2logical!(q::AbstractVector)
    q .= q[SA[1, 2, 7, 3, 4, 5, 6]]
end

function logical2abb!(q::AbstractVector)
    q .= q[SA[1, 2, 4, 5, 6, 7, 3]]
end

abb2logical(q::AbstractMatrix) = q[:,SA[1, 2, 7, 3, 4, 5, 6]]
logical2abb(q::AbstractMatrix) = q[:,SA[1, 2, 4, 5, 6, 7, 3]]
abb2logical(q::AbstractVector) = q[SA[1, 2, 7, 3, 4, 5, 6]]
logical2abb(q::AbstractVector) = q[SA[1, 2, 4, 5, 6, 7, 3]]
