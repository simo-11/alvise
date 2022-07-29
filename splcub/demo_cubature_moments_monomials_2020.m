
function demo_cubature_moments_monomials_2020

% This demo:
% * tests the performance of our new rules of degree ADE, on spline
%   curvilinear domains (only over random polynomials of total
%   degree at most ADE)
% * plots domain and pointsets
%
% Last update: June 4, 2020.

%============================= experiment data ============================

example=1;     % experiments: choose from 0 to 4 (see the attached file
%              define_example.m); for those in the paper
%              choose 1 or 2.

adeV=5:5:20;   % algebraic degree of precision

SPLtypestring='not-a-knot'; % cubic spline type (if necessary)

extraction_type=2; % Tchakaloff extraction:
% 1: lsqnonneg, 2: Fast Lawson-Hanson, 3: QR

number_experiments=1000; % number of tests for each degree.





%============================= main code below ============================

[polygon_sides,spline_order_vett]=define_example(example);

% ..................... reference routine ...................................
[nodes_x,nodes_y,weights]=splinegauss_2009b(max(adeV),polygon_sides,...
    [],[],[],spline_order_vett,0,SPLtypestring,4);

% ..................... new rule experiments ..............................

X=polygon_sides(:,1); Y=polygon_sides(:,2);

% Definition of the curvilinear sides via splines (Sx and SY)
[Sx,Sy]=Sx_Sy_preparation(X,Y,1,spline_order_vett,SPLtypestring);

for ii=1:length(adeV)
    
    % Rule ADE (algebraic degree of precision)
    n=adeV(ii);
    
    % Cubature nodes and weights stored in the N x 3 matrix "xyw"
    [xyw,momsres(ii),Z,Zin,cmom,bbox] = splcub(n,Sx,Sy,extraction_type);
    
    Vxy= cvand(n,xyw(:,1:2),bbox);
    moms_splcub=((xyw(:,3))'*Vxy)';
    
    Vnodes= cvand(n,[nodes_x nodes_y],bbox);
    moms_splgss=(weights'*Vnodes)';
    
    %format short e; [cmom-moms_splcub cmom-moms_splgss]
    momerr_splcub(ii)=norm(cmom-moms_splcub,2);
    momerr_splgss(ii)=norm(cmom-moms_splgss,2);
    
    
    jj=1;
    
    for mm=0:n
        for nn=0:n-mm
            
            % Random polynomial integrand of degree "n"
            k1=rand(1); k2=rand(1); c1=rand(1);
            f=@(x,y) x.^mm.*y.^nn;
            
            % ............. reference result ............
            I_exact(jj) = exact_integral(mm,nn,Sx,Sy);
            
            %         % ............. reference value by splinegauss_2009b ............
            fnodes=feval(f,nodes_x,nodes_y);
            I_2009(jj)=weights'*fnodes;
            
            % ................ reference value by splcub ......................
            fnodes=feval(f,xyw(:,1),xyw(:,2));
            I_2019(jj)=xyw(:,3)'*fnodes;
            
            % ...................... statistics ...............................
            
            % ... Absolute error ...
            aeV(jj)=abs(I_exact(jj)-I_2019(jj));
            
            % ... Absolute error ...
            aeV_2009(jj)=abs(I_exact(jj)-I_2009(jj));
            
            % deltaVRE(jj)=abs(I_2019(jj)-I_2009(jj));
            
            % ... Relative error ...
            if abs(I_exact) > 10^(-9) | abs(I_exact) < 10^(+12)
                reV(jj)=aeV(jj)/abs(I_exact(jj));
                reV_2009(jj)=aeV_2009(jj)/abs(I_exact(jj));
            else
                reV(jj)=0;
                reV_2009(jj)=0;
            end
            
            jj=jj+1;
        end
    end
    
    % absolute and relative errors
    aeVmax(ii)=max(aeV);
    [reVmax(ii),imaxR(ii)]=max(reV);
    
    aeVmax_2009(ii)=max(aeV_2009);
    reVmax_2009(ii)=max(reV_2009);
    
    % deltaVmaxRE(ii)=max(deltaVRE);
    
    % mean (log based)
    iok1=find(aeV > 0);
    log_aeV(ii)=10^(mean(log10(aeV(iok1))));
    
    iok2=find(reV > 0);
    log_reV(ii)=10^(mean(log10(reV(iok2))));
    
    % mean (log based)
    iok3=find(aeV_2009 > 0);
    log_aeV_2009(ii)=10^(mean(log10(aeV_2009(iok3))));
    
    iok4=find(reV_2009 > 0);
    log_reV_2009(ii)=10^(mean(log10(reV_2009(iok4))));
    
    % mean 
    mean_abs(ii)=mean(aeV);
    mean_rel(ii)=mean(reV);
    
    % mean
    mean_abs_2009(ii)=mean(aeV_2009);
    mean_rel_2009(ii)=mean(reV_2009);
    
    fprintf('\n \t DEG: %2.0f max(RE) 2020: %1.1e 2009:%1.1e',...
        n,reVmax(ii),reVmax_2009(ii))
    fprintf('\n \t AE momerr splcub: %1.1e momerrsplcub: %1.1e \n \n',...
        momerr_splcub(ii),momerr_splgss(ii));
end

% .........................  display statistics ...........................

% statistics: domains and function
fprintf('\n \n \t domain: %2.0f',example);
switch extraction_type
    case 1
        fprintf('\n \t extraction type: lsqnonneg');
    case 2
        fprintf('\n \t extraction type: fast Lawson-Hanson');
    otherwise
        fprintf('\n \t extraction type: QR');
end

% statistics: header
fprintf('\n \n');
fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n \t 2020');
fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n  |  ade |  max(ae) |  max(re) | max*(ae) |  max*(re) |  mean_abs | mean_rel |');
fprintf('\n  ----------------------------------------------------------------------------');

% statistics: experiments
for ii=1:length(adeV)
    fprintf('\n  | %3.0f  | %1.2e | %1.2e | %1.2e | %1.2e  | %1.2e  | %1.2e |',...
        adeV(ii), aeVmax(ii),reVmax(ii),log_aeV(ii),log_reV(ii),...
        mean_abs(ii),mean_rel(ii));
end
fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n  ade      : algebraic degree of precision');
fprintf('\n  max(ae)  : worst absolute error');
fprintf('\n  max(re)  : worst relative error');
fprintf('\n  max*(ae) : log average absolute error');
fprintf('\n  max*(re) : log average relative error');
fprintf('\n  mean_abs : average absolute error');
fprintf('\n  mean_rel : average relative error');
fprintf('\n  ----------------------------------------------------------------------------');


fprintf('\n \n');
fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n \t 2009');
fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n  |  ade |  max(ae) |  max(re) | max*(ae) |  max*(re) |  mean_abs | mean_rel |');
fprintf('\n  ----------------------------------------------------------------------------');

% statistics: experiments
for ii=1:length(adeV)
    fprintf('\n  | %3.0f  | %1.2e | %1.2e | %1.2e | %1.2e  | %1.2e  | %1.2e |',...
        adeV(ii), aeVmax_2009(ii),reVmax_2009(ii),log_aeV_2009(ii),log_reV_2009(ii),...
        mean_abs_2009(ii),mean_rel_2009(ii));
end

fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n  ade      : algebraic degree of precision');
fprintf('\n  max(ae)  : worst absolute error');
fprintf('\n  max(re)  : worst relative error');
fprintf('\n  max*(ae) : log average absolute error');
fprintf('\n  max*(re) : log average relative error');
fprintf('\n  mean_abs : average absolute error');
fprintf('\n  mean_rel : average relative error');
fprintf('\n  ----------------------------------------------------------------------------');



fprintf('\n \n');
fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n \t MOMENTS');
fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n  |  ade |  moms_spc |  moms_spg |');
fprintf('\n  ----------------------------------------------------------------------------');

for ii=1:length(adeV)
    fprintf('\n  | %3.0f  | %1.2e  | %1.2e  |',...
        adeV(ii), momerr_splcub(ii),momerr_splgss(ii));
end

fprintf('\n  ----------------------------------------------------------------------------');
fprintf('\n  ade      : algebraic degree of precision');
fprintf('\n  moms_spc : momerr splcub 2 norm');
fprintf('\n  moms_spg : momerr splgauss 2 norm');
fprintf('\n  ----------------------------------------------------------------------------');



fprintf('\n \n');

% format short e
fprintf('\n moments in absolute value: min: %1.1e max: %1.1e \n \n',...
    min(abs(I_exact)),max(abs(I_exact)) );







function [polygon_sides,spline_order_vett]=define_example(example)

%--------------------------------------------------------------------------
% POLYGON DEFINITION.
%--------------------------------------------------------------------------
%
% INPUT:
%
% example: PARAMETER THAT CHOOSES THE EXAMPLE.
%
% OUTPUT:
% polygon_sides: DEFINES THE COORDINATES OF THE VERTICES (REPEATING AT
%                THE END THE FIRST POINT). THE VARIABLE "polygon_type" IS
%                DEFINED ABOVE;
%                "boundary_pts" IS A VECTOR HAVING AS COMPONENTS THE
%                VERTICES OF THE COMPONENT, DESCRIBED COUNTERCLOCKWISE.
%                POLYGON IS DEFINED INSIDE THE SQUARE [0,1] x [0,1].
%
% spline_order_vett: SPLINE TO BE USED IN SUBINTERVALS, DEFINED AS
%                "spline order" AND INDEX OF THE LAST POINT IN SUBINTERVAL.
%                THUS [4 5; 2 12] USES A SPLINE OF ORDER 4 (CUBIC) IN THE
%                INTERVAL DEFINED BY THE FIRST 5 POINTS, AND OF ORDER 2
%                (LINEAR SPLINE) IN THE INTERVAL FROM THE FIFTH POINT TO
%                THE LAST (I.E. 12th POINT).
%
%--------------------------------------------------------------------------
% OBSERVE THAT THE FIRST VERTEX IS REPEATED (TO CLOSE THE POLYGON).
%--------------------------------------------------------------------------

switch example
    
    case 1
        polygon_sides=[-1.000000000000000e+00 0.000000000000000e+00
            -2.000000000000000e+00 -1.000000000000000e+00
            -1.500000000000000e+00 -2.000000000000000e+00
            0.000000000000000e+00 -1.600000000000000e+00
            0.000000000000000e+00 -1.000000000000000e+00
            -0.2 -0.5
            -0.38 -0.75
            -0.20 -0.94
            -0.57 -1.28
            -1.000000000000000e+00 0.000000000000000e+00];
        
        % polygon_sides=[polygon_sides(:,1)+1 polygon_sides(:,2)+1];
        % spline_order_vett=[2 5; 4 size(polygon_sides,1)];
        polygon_sides=[polygon_sides(:,1)+1 polygon_sides(:,2)+1];
        spline_order_vett=[2 5; 4 size(polygon_sides,1)];
        
    case 2
        polygon_sides=[   0.49804687500000   0.23832417582418
            0.52148437500000   0.29601648351648
            0.55507812500000   0.35233516483516
            0.58398437500000   0.40315934065934
            0.58242187500000   0.42513736263736
            0.56757812500000   0.43475274725275
            0.54414062500000   0.41826923076923
            0.52226562500000   0.38942307692308
            0.49570312500000   0.38392857142857
            0.48710937500000   0.41552197802198
            0.48945312500000   0.47184065934066
            0.49882812500000   0.54601648351648
            0.50195312500000   0.62019230769231
            0.49257812500000   0.64766483516484
            0.47773437500000   0.63667582417582
            0.46601562500000   0.55975274725275
            0.45039062500000   0.47733516483516
            0.45273437500000   0.57760989010989
            0.45117187500000   0.65728021978022
            0.43085937500000   0.67788461538462
            0.41757812500000   0.65178571428571
            0.41679687500000   0.57074175824176
            0.41445312500000   0.48008241758242
            0.40195312500000   0.55700549450549
            0.38632812500000   0.63667582417582
            0.37304687500000   0.64766483516484
            0.36210937500000   0.62568681318681
            0.36992187500000   0.53640109890110
            0.37695312500000   0.44848901098901
            0.34492187500000   0.50068681318681
            0.32070312500000   0.54876373626374
            0.30195312500000   0.55151098901099
            0.29882812500000   0.53090659340659
            0.32460937500000   0.45398351648352
            0.34570312500000   0.40041208791209
            0.34960937500000   0.32074175824176
            0.35976562500000   0.23832417582418
            0.49804687500000   0.23832417582418];
        
        polygon_sides=[polygon_sides(:,1)-0.4 polygon_sides(:,2)-0.4];
        
        spline_order_vett=[4 size(polygon_sides,1)];
        
    case 3
        polygon_sides=[-1.000000000000000e+00 0.000000000000000e+00
            -2.000000000000000e+00 -1.000000000000000e+00
            -1.500000000000000e+00 -2.000000000000000e+00
            0.000000000000000e+00 -1.600000000000000e+00
            0.000000000000000e+00 -1.000000000000000e+00
            -1.570731731182077e-02 -9.998766324816606e-01
            -3.141075907812830e-02 -9.995065603657316e-01
            -4.710645070964257e-02 -9.988898749619700e-01
            -6.279051952931321e-02 -9.980267284282716e-01
            -7.845909572784557e-02 -9.969173337331280e-01
            -9.410831331851485e-02 -9.955619646030800e-01
            -1.097343110910457e-01 -9.939609554551796e-01
            -1.253332335643046e-01 -9.921147013144778e-01
            -1.409012319375829e-01 -9.900236577165575e-01
            -1.564344650402310e-01 -9.876883405951377e-01
            -1.719291002794096e-01 -9.851093261547739e-01
            -1.873813145857246e-01 -9.822872507286887e-01
            -2.027872953565124e-01 -9.792228106217657e-01
            -2.181432413965424e-01 -9.759167619387474e-01
            -2.334453638559051e-01 -9.723699203976767e-01
            -2.486898871648553e-01 -9.685831611286310e-01
            -2.638730499653733e-01 -9.645574184577980e-01
            -2.789911060392296e-01 -9.602936856769430e-01
            -2.940403252323042e-01 -9.557930147983300e-01
            -3.090169943749476e-01 -9.510565162951535e-01
            -3.239174181981495e-01 -9.460853588275453e-01
            -3.387379202452914e-01 -9.408807689542255e-01
            -3.534748437792570e-01 -9.354440308298674e-01
            -3.681245526846778e-01 -9.297764858882515e-01
            -3.826834323650903e-01 -9.238795325112865e-01
            -3.971478906347811e-01 -9.177546256839809e-01
            -4.115143586051092e-01 -9.114032766354451e-01
            -4.257792915650729e-01 -9.048270524660194e-01
            -4.399391698559154e-01 -8.980275757606155e-01
            -4.539904997395469e-01 -8.910065241883678e-01
            -4.679298142605735e-01 -8.837656300886934e-01
            -4.817536741017153e-01 -8.763066800438636e-01
            -4.954586684324074e-01 -8.686315144381913e-01
            -5.090414157503719e-01 -8.607420270039433e-01
            -5.224985647159486e-01 -8.526401643540924e-01
            -5.358267949789971e-01 -8.443279255020149e-01
            -5.490228179981321e-01 -8.358073613682701e-01
            -5.620833778521309e-01 -8.270805742745616e-01
            -5.750052520432788e-01 -8.181497174250233e-01
            -5.877852522924732e-01 -8.090169943749473e-01
            -6.004202253258841e-01 -7.996846584870905e-01
            -6.129070536529765e-01 -7.901550123756904e-01
            -6.252426563357051e-01 -7.804304073383298e-01
            -6.374239897486895e-01 -7.705132427757894e-01
            -6.494480483301841e-01 -7.604059656000306e-01
            -6.613118653236519e-01 -7.501110696304595e-01
            -6.730125135097736e-01 -7.396310949786095e-01
            -6.845471059286889e-01 -7.289686274214113e-01
            -6.959127965923145e-01 -7.181262977631887e-01
            -7.071067811865477e-01 -7.071067811865475e-01
            -7.181262977631889e-01 -6.959127965923143e-01
            -7.289686274214116e-01 -6.845471059286887e-01
            -7.396310949786099e-01 -6.730125135097731e-01
            -7.501110696304597e-01 -6.613118653236517e-01
            -7.604059656000310e-01 -6.494480483301835e-01
            -7.705132427757893e-01 -6.374239897486896e-01
            -7.804304073383297e-01 -6.252426563357052e-01
            -7.901550123756906e-01 -6.129070536529763e-01
            -7.996846584870907e-01 -6.004202253258839e-01
            -8.090169943749475e-01 -5.877852522924730e-01
            -8.181497174250235e-01 -5.750052520432785e-01
            -8.270805742745618e-01 -5.620833778521306e-01
            -8.358073613682704e-01 -5.490228179981315e-01
            -8.443279255020152e-01 -5.358267949789964e-01
            -8.526401643540923e-01 -5.224985647159487e-01
            -8.607420270039439e-01 -5.090414157503709e-01
            -8.686315144381914e-01 -4.954586684324072e-01
            -8.763066800438637e-01 -4.817536741017150e-01
            -8.837656300886936e-01 -4.679298142605732e-01
            -8.910065241883679e-01 -4.539904997395467e-01
            -8.980275757606156e-01 -4.399391698559151e-01
            -9.048270524660195e-01 -4.257792915650727e-01
            -9.114032766354452e-01 -4.115143586051089e-01
            -9.177546256839813e-01 -3.971478906347804e-01
            -9.238795325112868e-01 -3.826834323650897e-01
            -9.297764858882515e-01 -3.681245526846779e-01
            -9.354440308298675e-01 -3.534748437792568e-01
            -9.408807689542256e-01 -3.387379202452911e-01
            -9.460853588275454e-01 -3.239174181981492e-01
            -9.510565162951536e-01 -3.090169943749473e-01
            -9.557930147983301e-01 -2.940403252323039e-01
            -9.602936856769431e-01 -2.789911060392293e-01
            -9.645574184577982e-01 -2.638730499653726e-01
            -9.685831611286312e-01 -2.486898871648546e-01
            -9.723699203976767e-01 -2.334453638559053e-01
            -9.759167619387474e-01 -2.181432413965425e-01
            -9.792228106217657e-01 -2.027872953565125e-01
            -9.822872507286887e-01 -1.873813145857243e-01
            -9.851093261547740e-01 -1.719291002794093e-01
            -9.876883405951378e-01 -1.564344650402307e-01
            -9.900236577165575e-01 -1.409012319375826e-01
            -9.921147013144779e-01 -1.253332335643043e-01
            -9.939609554551797e-01 -1.097343110910454e-01
            -9.955619646030800e-01 -9.410831331851410e-02
            -9.969173337331280e-01 -7.845909572784482e-02
            -9.980267284282716e-01 -6.279051952931335e-02
            -9.988898749619700e-01 -4.710645070964227e-02
            -9.995065603657316e-01 -3.141075907812799e-02
            -9.998766324816606e-01 -1.570731731182046e-02
            -1.000000000000000e+00 0.000000000000000e+00];
        
        spline_order_vett=[2 5; 4 size(polygon_sides,1)];
        
    case 4
        polygon_sides=[-1.000000000000000e+00 0.000000000000000e+00
            -2.000000000000000e+00 -1.000000000000000e+00
            -1.500000000000000e+00 -2.000000000000000e+00
            0.000000000000000e+00 -1.600000000000000e+00
            0.000000000000000e+00 -1.000000000000000e+00
            0.000000000000000e+00 -1.000000000000000e+00
            -1.233675183394123e-04 -9.842926826881794e-01
            -4.934396342684000e-04 -9.685892409218717e-01
            -1.110125038029985e-03 -9.528935492903573e-01
            -1.973271571728441e-03 -9.372094804706866e-01
            -3.082666266872036e-03 -9.215409042721551e-01
            -4.438035396920004e-03 -9.058916866814857e-01
            -6.039044544820293e-03 -8.902656889089547e-01
            -7.885298685522124e-03 -8.746667664356957e-01
            -9.976342283442463e-03 -8.590987680624174e-01
            -1.231165940486223e-02 -8.435655349597692e-01
            -1.489067384522613e-02 -8.280708997205904e-01
            -1.771274927131128e-02 -8.126186854142754e-01
            -2.077718937823425e-02 -7.972127046434875e-01
            -2.408323806125257e-02 -7.818567586034575e-01
            -2.763007960232344e-02 -7.665546361440947e-01
            -3.141683887136892e-02 -7.513101128351451e-01
            -3.544258154220192e-02 -7.361269500346270e-01
            -3.970631432305693e-02 -7.210088939607707e-01
            -4.420698520166988e-02 -7.059596747676961e-01
            -4.894348370484647e-02 -6.909830056250525e-01
            -5.391464117245470e-02 -6.760825818018505e-01
            -5.911923104577455e-02 -6.612620797547086e-01
            -6.455596917013262e-02 -6.465251562207429e-01
            -7.022351411174854e-02 -6.318754473153221e-01
            -7.612046748871326e-02 -6.173165676349102e-01
            -8.224537431601886e-02 -6.028521093652194e-01
            -8.859672336455471e-02 -5.884856413948912e-01
            -9.517294753398042e-02 -5.742207084349273e-01
            -1.019724242393844e-01 -5.600608301440848e-01
            -1.089934758116321e-01 -5.460095002604533e-01
            -1.162343699113065e-01 -5.320701857394267e-01
            -1.236933199561364e-01 -5.182463258982847e-01
            -1.313684855618088e-01 -5.045413315675924e-01
            -1.392579729960564e-01 -4.909585842496287e-01
            -1.473598356459077e-01 -4.775014352840512e-01
            -1.556720744979849e-01 -4.641732050210033e-01
            -1.641926386317297e-01 -4.509771820018682e-01
            -1.729194257254382e-01 -4.379166221478694e-01
            -1.818502825749766e-01 -4.249947479567214e-01
            -1.909830056250525e-01 -4.122147477075269e-01
            -2.003153415129094e-01 -3.995797746741161e-01
            -2.098449876243096e-01 -3.870929463470235e-01
            -2.195695926616702e-01 -3.747573436642949e-01
            -2.294867572242107e-01 -3.625760102513104e-01
            -2.395940343999691e-01 -3.505519516698163e-01
            -2.498889303695404e-01 -3.386881346763482e-01
            -2.603689050213903e-01 -3.269874864902267e-01
            -2.710313725785884e-01 -3.154528940713114e-01
            -2.818737022368111e-01 -3.040872034076857e-01
            -2.928932188134524e-01 -2.928932188134525e-01
            -3.040872034076857e-01 -2.818737022368112e-01
            -3.154528940713113e-01 -2.710313725785884e-01
            -3.269874864902267e-01 -2.603689050213903e-01
            -3.386881346763481e-01 -2.498889303695405e-01
            -3.505519516698163e-01 -2.395940343999691e-01
            -3.625760102513103e-01 -2.294867572242107e-01
            -3.747573436642947e-01 -2.195695926616703e-01
            -3.870929463470235e-01 -2.098449876243096e-01
            -3.995797746741159e-01 -2.003153415129095e-01
            -4.122147477075268e-01 -1.909830056250527e-01
            -4.249947479567214e-01 -1.818502825749766e-01
            -4.379166221478693e-01 -1.729194257254382e-01
            -4.509771820018683e-01 -1.641926386317297e-01
            -4.641732050210035e-01 -1.556720744979849e-01
            -4.775014352840511e-01 -1.473598356459078e-01
            -4.909585842496288e-01 -1.392579729960564e-01
            -5.045413315675924e-01 -1.313684855618087e-01
            -5.182463258982848e-01 -1.236933199561363e-01
            -5.320701857394267e-01 -1.162343699113065e-01
            -5.460095002604531e-01 -1.089934758116322e-01
            -5.600608301440849e-01 -1.019724242393844e-01
            -5.742207084349273e-01 -9.517294753398042e-02
            -5.884856413948911e-01 -8.859672336455482e-02
            -6.028521093652195e-01 -8.224537431601886e-02
            -6.173165676349102e-01 -7.612046748871326e-02
            -6.318754473153219e-01 -7.022351411174865e-02
            -6.465251562207428e-01 -6.455596917013273e-02
            -6.612620797547085e-01 -5.911923104577455e-02
            -6.760825818018505e-01 -5.391464117245470e-02
            -6.909830056250525e-01 -4.894348370484647e-02
            -7.059596747676959e-01 -4.420698520166988e-02
            -7.210088939607705e-01 -3.970631432305705e-02
            -7.361269500346272e-01 -3.544258154220192e-02
            -7.513101128351453e-01 -3.141683887136892e-02
            -7.665546361440945e-01 -2.763007960232344e-02
            -7.818567586034573e-01 -2.408323806125268e-02
            -7.972127046434873e-01 -2.077718937823425e-02
            -8.126186854142753e-01 -1.771274927131139e-02
            -8.280708997205904e-01 -1.489067384522613e-02
            -8.435655349597690e-01 -1.231165940486223e-02
            -8.590987680624171e-01 -9.976342283442463e-03
            -8.746667664356955e-01 -7.885298685522235e-03
            -8.902656889089546e-01 -6.039044544820293e-03
            -9.058916866814857e-01 -4.438035396920004e-03
            -9.215409042721550e-01 -3.082666266872036e-03
            -9.372094804706865e-01 -1.973271571728441e-03
            -9.528935492903573e-01 -1.110125038029985e-03
            -9.685892409218716e-01 -4.934396342684000e-04
            -9.842926826881794e-01 -1.233675183394123e-04
            -1.000000000000000e+00 0.000000000000000e+00];
        
        spline_order_vett=[2 5; 4 size(polygon_sides,1)];
        
    case 5
        polygon_sides=[-1.000000000000000e+00 -1.000000000000000e+00
            +1.000000000000000e+00 -1.000000000000000e+00
            +1.000000000000000e+00 +1.000000000000000e+00
            -1.000000000000000e+00 +1.000000000000000e+00
            -1.000000000000000e+00 -1.000000000000000e+00];
        
        spline_order_vett=[2 size(polygon_sides,1)];
        
        
        
end



function [f,fstr]=cubature_gallery(function_type,parms)

%--------------------------------------------------------------------------
% FUNCTION DEFINITION.
%--------------------------------------------------------------------------
%
% INPUT:
%
% function_type: PARAMETER THAT CHOOSES THE FUNCTION.
% parms: SPECIFIC VALUES TO CHANGE THE TEST FUNCTION.
%
% OUTPUT:
%
% f: function to test
% fstr: string representing the function
%--------------------------------------------------------------------------

if nargin < 2,
    k=1; x0=0; y0=0;
else
    L=length(parms);
    
    switch L
        case 0
            k=1; x0=0; y0=0;
        case 1
            k=parms(1); x0=0; y0=0;
        case 2
            k=1; x0=parms(1); y0=parms(2);
        case 3
            k=parms(1); x0=parms(2); y0=parms(3);
    end
    
end

ks=num2str(k); x0s=num2str(x0); y0s=num2str(y0);


switch function_type
    
    case 1
        fstr=strcat('exp(-',ks,'*((x-',x0s,').^2+(y-',y0s,').^2))');
        f=@(x,y) exp(-k*((x-x0).^2+(y-y0).^2));
    case 2 % low regular, e.g. set k=1/2, k=11/2
        fstr=strcat('( (x-',x0s,').^2 +(y-',y0s,').^2 ).^',ks);
        f=@(x,y)( (x-x0).^2 +(y-y0).^2 ).^k;
    case 3
        fstr=strcat('cos(',ks,'*(x+y))');
        f=@(x,y) cos(k*(x+y));
    case 4 % fixed polynomial
        fstr=strcat('((x-',x0s,')+(y-',y0s,')).^',ks);
        f=@(x,y) ((x-x0)+(y-y0)).^k;
    case 5
        fstr=strcat('exp(-',ks,'*((x-',x0s,').^2+(y-',y0s,').^2))');
        f=@(x,y) exp(-k*((x-x0).^2+(y-y0).^2));
    case 6 % Franke
        fstr='franke';
        f=@(x,y) .75*exp(-((9*x-2).^2 + (9*y-2).^2)/4) + ...
            .75*exp(-((9*x+1).^2)/49 - (9*y+1)/10) + ...
            .5*exp(-((9*x-7).^2 + (9*y-3).^2)/4) - ...
            .2*exp(-(9*x-4).^2 - (9*y-7).^2);
    case 7 % exponential
        fstr=strcat('exp(',k,'*(x-y))');
        f=@(x,y) exp(k*(x-y));
    case 8
        fstr='ones(size(x))';
        f=@(x,y) ones(size(x));
    case 9 % polynomial of degree "degree"
        fstr='nonstandard polynomial';
        f=@(x,y) (1+7.671018281063714e-01*(x-x0) + ...
            3.547875152618030e-02*(y-y0)).^k;
    case 10 % polynomial of degree "degree"
        fstr='random polynomial';
        f=@(x,y) (rand(1)+rand(1)*(x-x0) + ...
            rand(1)*(y-y0)).^k;
end




function [Sx,Sy,Nsub]=Sx_Sy_preparation(X,Y,Nsub,spline_parms,SPLtypestring)

% OBJECT:
%
%
% INPUT:
% X,Y:
% Nsub:
% spline_parms:
% SPLtypestring:
%
% OUTPUT:
% Sx,Sy:
% Nsub:
%

% check input variables.
L=size(spline_parms,1);

for i=1:L % define spline
    
    if i == 1
        imin=1;
    else
        imin=spline_parms(i-1,2);
    end
    imax=spline_parms(i,2);
    
    tL=imin:imax;
    xL=X(imin:imax);
    yL=Y(imin:imax);
    
    [SxL,SyL]=compute_parametric_spline(tL,xL,yL,...
        spline_parms(i,1),SPLtypestring);
    
    Sx(i)=SxL;
    Sy(i)=SyL;
    
end





%--------------------------------------------------------------------------
% compute_parametric_spline
%--------------------------------------------------------------------------

function [ppx,ppy]=compute_parametric_spline(s,x,y,spline_order,...
    SPLtypestring)

% OBJECT:
% compute parametric spline relavant parameters "ppx", "ppy" so that a
% point at the boundary of the  domain has coordinates (x(s),y(s))
%
% INPUT:
% s: parameter data.
% x: determine spline x(s) interpolating (s,x)
% y: determine spline y(s) interpolating (s,y)
% spline_order: spline order (i.e. degree + 1)
% SPLtypestring: string with the spline type i.e.
%             'complete'   : match endslopes (as given in VALCONDS, with
%                     default as under *default*).
%             'not-a-knot' : make spline C^3 across first and last interior
%                     break (ignoring VALCONDS if given).
%             'periodic'   : match first and second derivatives at first
%                     data point with those at last data point (ignoring
%                     VALCONDS if given).
%             'second'     : match end second derivatives (as given in
%                    VALCONDS, with default [0 0], i.e., as in variational).
%             'variational': set end second derivatives equal to zero
%                     (ignoring VALCONDS if given).
%
% OUTPUT:
% ppx: spline x(t) data
% ppy: spline y(t) data


switch spline_order
    case 4
        
        % CUBIC SPLINES BY CSAPE. AS DEFAULT WE USE PERIODIC CUBIC SPLINES.
        % Derivatives parameters are computed as well.
        ppx=csape(s,x,SPLtypestring);
        ppy=csape(s,y,SPLtypestring);
        
    otherwise
        
        ppx=spapi(spline_order,s,x);
        ppy=spapi(spline_order,s,y);
        
end

if (ppx.form =='B-')
    ppx=fn2fm(ppx,'pp');
    ppy=fn2fm(ppy,'pp');
end





%--------------------------------------------------------------------------
% Attached functions.
%--------------------------------------------------------------------------

function plot_spl(Sx,Sy,P,Z,Zin)

% ........................ plot preferences ...............................
plot_axis=1;

% ........................... plot: start .................................
clf;
hold on;

if plot_axis == 0, axis off; else, axis on; end

% ........................ plot: grids .........................

plot(Z(:,1),Z(:,2),'ro');
plot(Zin(:,1),Zin(:,2),'go');
plot(P(:,1),P(:,2),'ko','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',2)

% ....................... plot: parametric splines ........................

L=length(Sx);
for ii=1:L
    
    SxL=Sx(ii);
    SyL=Sy(ii);
    
    SxL_breaks=SxL.breaks;
    
    
    Nbreaks=length(SxL_breaks);
    
    for kk=1:Nbreaks-1
        
        t0=SxL_breaks(kk); t1=SxL_breaks(kk+1);
        
        N=10^3;
        tt=linspace(t0,t1,N);
        
        xx=ppval(SxL,tt);
        yy=ppval(SyL,tt);
        
        
        plot(xx,yy,'k-','LineWidth',1);
    end
    
end

hold off;






%==========================================================================
% cvand
%==========================================================================

function V = cvand(n,pts,bbox)

% OBJECT:
% computes by recurrence the Chebyshev-Vandermonde matrix on a 2d
% arbitrarily located mesh, in a total-degree product Chebyshev basis
% of a given rectangle
%
% INPUT:
% n: polynomial degree
% pts: 2-column array of mesh point coordinates
% bbox: [bbox(1),bbox(2)] x [bbox(3),bbox(4)] bounding box for the mesh
%
% OUTPUT:
% V: Chebyshev-Vandermonde matrix (V is a matrix "M x N" where "M" is the
% number of points and "N" the dimension of the polynomial space, i.e.
%                     N=(M+1)*(M+2)/2
%
% DATA:
%
% built: March 2019
% check: June 4, 2020


% default rectangle containing the mesh if not passed
if isempty(bbox)
    bbox=[min(pts(:,1)) max(pts(:,1)) min(pts(:,2)) max(pts(:,2))];
end

a=bbox(1);b=bbox(2);c=bbox(3);d=bbox(4);

Vx=chebpolys(n,(2*pts(:,1)-b-a)/(b-a));
Vy=chebpolys(n,(2*pts(:,2)-d-c)/(d-c));

k=0;
for i=0:n
    for j=0:n-i
        k=k+1;
        V(:,k)=Vx(:,i+1).*Vy(:,j+1);
    end
end




%==========================================================================
% chebpolys
%==========================================================================

function T=chebpolys(deg,x)

% OBJECT:
% computes the Chebyshev-Vandermonde matrix on the real line by recurrence
%
% INPUT:
% deg = maximum polynomial degree
% x = 1-column array of abscissas
%
% OUTPUT
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg
%
% DATA:
%
% built: March 2019
% check: May 2, 2019

% inizialization
T=zeros(length(x),deg+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

% 3-term recurrence
for j=2:deg
    t2=2*x.*t1-t0;
    T(:,j+1)=t2;
    t0=t1;
    t1=t2;
end





























%==========================================================================
% chebmom
%==========================================================================

function [Iexact] = exact_integral(mm,nn,Sx,Sy)

% OBJECT:
% computes the exact integral up to degree n of a total-degree polynomial
% x.^mm.*y.^nn to the Lebesgue measure in a Jordan spline polygon,
% whose boundary is given by the counterclockwise concatened spline
% arcs (Sx,Sy)
%
% INPUT:
% mm,nn: constants for polynomial definition
% Sx,Sy: arrays of spline structures; (SX(i),Sy(i)) is the i-th arc
% the arcs must be counterclockwise concatenated forming a Jordan curve
% Sx(i).breaks(end)=Sx(i+1).breaks(1),Sy(i).breaks(end)=Sy(i+1).breaks(1)
% i=1,...,end, with end+1:=1

% OUTPUT:
% Iexact: integral value
%
% DATA:
% built: June 6, 2020
% check: June 6, 2020
% modified: June 6, 2020

n=mm+nn;
xyw=lineint(n,Sx,Sy);
intV=intpoly(mm,nn,xyw(:,1:2));
Iexact=intV'*xyw(:,3);




%==========================================================================
% lineint
%==========================================================================

function xyw = lineint(m,Sx,Sy)

% OBJECT:
% computes weights and nodes for the line integral on concatenated
% spline arcs; the formula is exact on bivariate polynomials up to deg m
%
% INPUT:
% m = polynomial degree
% Sx,Sy: arrays of spline structures; (SX(i),Sy(i)) is the i-th arc
% the arcs must be concatenated
% Sx(i).breaks(end)=Sx(i+1).breaks(1),Sy(i).breaks(end)=Sy(i+1).breaks(1)
%
% OUTPUT:
% xyw: 3-column array of nodes coords (cols 1-2) and weights (col 3)
%
% DATA:
% built: March 2019
% check: May 2, 2019

xyw=[];
for i=1:length(Sx)
    
    ord_spline=Sx(i).order; % degree is the spline order minus 1.
    
    % Gauss-Legendre nodes in [-1,1] and corresponding weights on the side.
    k=ceil(((ord_spline-1)*(m+2))/2);
    ab=r_jacobi(k,0,0);
    xw=gauss(k,ab);
    t=xw(:,1); w=xw(:,2);
    
    a=Sx(i).breaks(1:end-1); b=Sx(i).breaks(2:end);
    alpha=(b-a)/2; beta=(b+a)/2;
    dSy(i)=fnder(Sy(i));
    for j=1:length(a)
        nodes=alpha(j)*t+beta(j);
        wloc=w*alpha(j);
        xyw=[xyw;[ppval(Sx(i),nodes) ppval(Sy(i),nodes) ...
            wloc.*ppval(dSy(i),nodes)]];
    end
end





%==========================================================================
% intcvand
%==========================================================================

function intV = intpoly(mm,nn,pts)

% computes by recurrence an x-primitive of the polynomial x^mm * y^nn on a
% 2d arbitrarily located mesh

% June 2020

% INPUT:
% mm,nn: polynomial parameters
% pts: 2-column array of mesh point coordinates

% OUTPUT:
% intV: x-primitive of (a+bx+cy)^n evaluation at pts (column vector)

x=pts(:,1);
y=pts(:,2);

intV=(1/(mm+1))*(x.^(mm+1)).*y.^(nn);












%==========================================================================
% r_jacobi
%==========================================================================

function ab=r_jacobi(N,a,b)
%R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%   AB=R_JACOBI(N,A,B) generates the Nx2 array AB of the first
%   N recurrence coefficients for the monic Jacobi polynomials
%   orthogonal on [-1,1] relative to the weight function
%   w(x)=(1-x)^A*(1+x)^B. The call AB=R_JACOBI(N,A) is the same
%   as AB=R_JACOBI(N,A,A) and AB=R_JACOBI(N) the same as
%   AB=R_JACOBI(N,0,0).
%
%   Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%   Gautschi, 4-4-2002.
%   Edited by Walter Leopardi 10-22-2006.

if nargin<2, a=0; end
if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
nu=(b-a)/(a+b+2);
if a+b+2 > 128
    mu=exp((a+b+1)*log(2)+((gammaln(a+1)+gammaln(b+1))-gammaln(a+b+2)));
else
    mu=2^(a+b+1)*((gamma(a+1)*gamma(b+1))/gamma(a+b+2));
end
if N==1, ab=[nu mu]; return, end
N=N-1; n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];





%==========================================================================
% gauss
%==========================================================================

function xw=gauss(N,ab)
%GAUSS Gauss quadrature rule.
%   GAUSS(N,AB) generates the Nx2 array XW of Gauss quadrature
%   nodes and weights for a given weight function W. The nodes,
%   in increasing order, are placed into the first column of XW,
%   and the corresponding weights into the second column. The
%   weight function W is specified by the Nx2 input array AB
%   of recurrence coefficients for the polynomials orthogonal
%   with respect to the weight function W.

N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
    J(n,n-1)=sqrt(ab(n,2));
    J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];







