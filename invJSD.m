function [JSD]=invJSD(P,Q)
% gives 1-JSD (similarity, not distance) 

% make sure vectors are normalized (values add up to 1)
if sum(P)~=1
    P=P./sum(P);
end

if sum(Q)~=1
    Q=Q./sum(Q);
end

M=0.5*(P + Q);

P_M=(P.*log(P./M));

P_M(isnan(P_M))=[];
P_M(isinf(P_M))=[];

Q_M=(Q.*log(Q./M));

Q_M(isnan(Q_M))=[];
Q_M(isinf(Q_M))=[];

JSD=1-sqrt(0.5 *sum(P_M) + 0.5 *sum(Q_M));
JSD(JSD==1)=nan;

end