%______________________________________________________________________________________________
%   (MVDE)
%
%  Author and programmer: Al-Attar Ali Mohamed
%
%         e-Mail: engatar@yahoo.com
%_______________________________________________________________________________________________
% The initial parameters that you need are:
clear all

% This algorithm is eligible to solve unimodal and multimodal optimization problems using a small population size
for iji=1:23
    
    if iji==1;F=('F1');elseif iji==2;F=('F2');elseif iji==3;F=('F3');elseif iji==4;F=('F4');elseif iji==5;F=('F5'); ...
    elseif iji==6;F=('F6');elseif iji==7; F=('F7'); elseif iji==8; F=('F8');elseif iji==9; F=('F9'); ...
    elseif iji==10; F=('F10');elseif iji==11; F=('F11');elseif iji==12; F=('F12'); ...
    elseif iji==13; F=('F13');elseif iji==14; F=('F14');elseif iji==15; F=('F15');
    elseif iji==16; F=('F16');elseif iji==17; F=('F17');elseif iji==18; F=('F18');
    elseif iji==19; F=('F19');elseif iji==20; F=('F20');elseif iji==21; F=('F21');
    elseif iji==22; F=('F22');elseif iji==23; F=('F23');
    end
    
    if iji < 14;MAX_ITER=1000;else; MAX_ITER=500;end% Maximum number of iterations
    
    n =20;            % Number of search agents
    s1=[];
for j=1:7    
    % Load details of the selected benchmark function
    [lb,ub,d,fobj] = Get_Functions_details(F);
    
    [Best_pos,Best_score,Convergence_curve]=MPSOGA(n,MAX_ITER,ub,lb,d,fobj);
    s1(j)=Best_score;

    %Draw and display objective function
    
    % figure,semilogy(Convergence_curve); title( F ); xlabel('Iteration'); ylabel('Best score obtained so far');
end

ss1=[min(s1),mean(s1), max(s1), std(s1)];
    display(['      min ','      mean ','     max ','     std ']);
    display(['SS_MPSOGA ',num2str(ss1)]);

end



