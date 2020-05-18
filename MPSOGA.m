%_________________________________________________________________________
% MPSOGA           
% 
%  Developed in MATLAB R2017b 
% 
%  Author and programmer: Al-Attar Ali Mohamed
%                                      
%         e-Mail: engatar@yahoo.com  
%                 attar@aswu.edu.eg 
%_________________________________________________________________________
% 

function [best,fmin,Convergence_curve]= MPSOGA(n,N_iter,Ub,Lb,d,Fun)

if length(Lb) < d 
    Lb=Lb*ones(1,d);Ub=Ub*ones(1,d);
end

% Initialize the population/solutions
    for i=1:d
        Sol(:,i)=rand(n,1).*(Ub(i)-Lb(i))+Lb(i);
    end
    for i=1:n
        Fitness(i)=Fun(Sol(i,:));  
        %     nfun=nfun+1;
    end
[Fitness, location]=sort(Fitness);
Sol=Sol(location,:);  
% Find the global best solution
fmin = Fitness(1);    best = Sol(1,:); 
% display(['At iteration: ', num2str(1), ' Best fitness is: ',num2str(fmin)]);

Convergence_curve(1)=fmin;
% Start the iterations  
for t=2:N_iter

%________________________________(  GA  )____________________________
keep=floor(0.5*n); % 
    [NewSol]=GA(n,Sol,Lb,Ub,d,keep);NewSol=NewSol(1:n,:); %
    NewSol=max(NewSol,repmat(Lb,n,1));
    NewSol=min(NewSol,repmat(Ub,n,1));    

    % Evaluate new solutions
    for i=keep+1:n
        Fnew(i)=Fun(NewSol(i,:));  
        %     nfun=nfun+1;
    end
%________________________________( MSPO  )____________________________
for i=1:keep%#ok<*ALIGN> 
    I=randperm(round(n/4));I(I==1)=[];   I=I(1);  %mating vectors  
%     I=randperm(round(n/3));I(I==1)=[];   I=I(1);  %mating vectors  


%     NewSol(i,:)=Sol(i,:)+0.01*(1-t/N_iter)*unifrnd(Lb-Sol(i,:),Ub-Sol(i,:))...
%         +(2*t/N_iter)*rand(size(best)).*(best-Sol(i,:))...
%         +(1-t/N_iter)*rand(size(best)).*(Sol(I,:)-Sol(i,:));%step away like PSO

    NewSol(i,:)=Sol(i,:)+0.05*(1-t/N_iter)*unifrnd(Lb-Sol(i,:),Ub-Sol(i,:))...
        +0.5*(2*t/N_iter)*rand(size(best)).*(best-Sol(i,:))...
        +2*(1-t/N_iter)*rand(size(best)).*(Sol(I,:)-Sol(i,:));%step away like PSO
   
    NewSol(i,:)=max(NewSol(i,:),Lb);
    NewSol(i,:)=min(NewSol(i,:),Ub);    
end    
    % Evaluate new solutions
    for i=1:keep
        Fnew(i)=Fun(NewSol(i,:));  
        %     nfun=nfun+1;
    end    
    % select the best solution
    [Fitness,I]=unique([Fitness, Fnew],'first'); 
    Fitness=Fitness(1:n);
    
    Sol=[Sol; NewSol]; 
    Sol=Sol(I(1:n),:);
    
    % Find the global best solution
    fmin = Fitness(1);    best = Sol(1,:);    
    Convergence_curve(t)=fmin;
    
%     if mod(t,1)==0 
%         disp(['At iteration: ', num2str(t), ' Best fitness is: ',num2str(fmin)]);
%     end
        
end



    function [par]=GA(popsize,par,lb,ub,dim,keep)
%_______________________________________________________
% III. GA parameters
mutrate=0.05*1;
nmut=ceil((popsize-1)*dim*mutrate); % total number of mutations
if size(ub,2)==1 % numnber of boundaries
    lb=lb*ones(1,dim);ub=ub*ones(1,dim);
end
%_______________________________________________________
% Pair and mate
M=ceil((popsize-keep)/2); % number of matings
prob=flipud([1:keep]'/sum([1:keep])); %#ok<NBRAK> % weights chromosomes
odds=[0 cumsum(prob(1:keep))']; % probability distribution function
pick1=rand(1,M); % mate #1
pick2=rand(1,M); % mate #2
% ma and pa contain the indicies of the chromosomes that will mate
ic=1;
while ic<=M;for id=2:keep+1;if pick1(ic)<=odds(id) && pick1(ic)>odds(id-1);ma(ic)=id-1; end;if pick2(ic)<=odds(id) && pick2(ic)>odds(id-1);pa(ic)=id-1; end;end;ic=ic+1;end %#ok<*AGROW>
%_______________________________________________________
% Performs mating using single point crossover
ix=1:2:keep; % index of mate #1
xp=ceil(rand(1,M)*dim); % crossover point
r=rand(1,M); % mixing parameter
for ic=1:M
xy=par(ma(ic),xp(ic))-par(pa(ic),xp(ic)); % ma and pa mate
par(keep+ix(ic),:)=par(ma(ic),:); % 1st offspring
par(keep+ix(ic)+1,:)=par(pa(ic),:); % 2nd offspring
par(keep+ix(ic),xp(ic))=par(ma(ic),xp(ic))-r(ic).*xy;% 1st
par(keep+ix(ic)+1,xp(ic))=par(pa(ic),xp(ic))+r(ic).*xy;% 2nd
if xp(ic)<dim % crossover when last variable not selected
par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic)) par(keep+ix(ic)+1,xp(ic)+1:dim)];
par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic)) par(keep+ix(ic),xp(ic)+1:dim)];
end % if
end
%_______________________________________________________
% Mutate the population
mrow=sort(ceil(rand(1,nmut)*(popsize/2))+popsize/2);mcol=ceil(rand(1,nmut)*dim);
for ii=1:nmut
par(mrow(ii),mcol(ii))=(ub(mcol(ii))-lb(mcol(ii)))*rand+lb(mcol(ii));% mutation
end% ii%__________________ End GA    ________________________

