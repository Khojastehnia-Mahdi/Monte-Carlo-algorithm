% Monte-Carlo algorithm
% capacity under the joint TPC+PAC constraints for a MIMO channel
% inputs: 
% W : Channel Gram matrix
% PT : total transmit power constraint
% P1 : per-antenna power constraints


function MC_TPC_PAC(W,PT,P1)

[Nt,~]=size(W);

% initial Tx covariance matrix
for i=1:Nt
    initial_entry_R(i)=min(P1(i),PT/Nt);
end
R{1}=diag(initial_entry_R);
R_star=R{1};
capacity(1)=log(det(eye(Nt,Nt)+W*R{1}));

% number of trials
number_of_trials=1e3;

for j=2:number_of_trials
	clear H
	var=1;
    
    % H is a random matrix
	H=sqrt(var)*(randn(Nt,Nt));
	bb=(H'*H);
    
    % random feasible transmit covariance
	for i=1:Nt
        rii(i)=bb(i,i);
        a(i)=(P1(i)/(rii(i)));
    end
	a=min(a);
	aa=PT/trace(H'*H);
    a=min(a,aa);
    R{j}=a.*H'*H;
    
    % The best transmit covariance
	C(j)=log(det(eye(Nt,Nt)+W*R{j}));
	if C(j)>capacity(j-1)
        R_star=R{j};
	end
        
    capacity(j)=log(det(eye(Nt,Nt)+W*R_star));
end
save('PAC_TPC_MC.mat')
end

