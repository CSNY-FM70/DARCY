function [uh,A_int,rhs]=twod_linear_FEM(ns,xv,yv,elt2vert,nvtx,ne,h,a,f)
    [Jks,invJks,detJks]=get_jac_info(xv,yv,ne,elt2vert);
    [Aks,bks]=get_elt_arrays2D(xv,yv,invJks,detJks,ne,elt2vert,a,f);
    A = sparse(nvtx,nvtx); b = zeros(nvtx,1);
    for row_no=1:3
        nrow=elt2vert(:,row_no);
        for col_no=1:3
            ncol=elt2vert(:,col_no);
            A=A+sparse(nrow,ncol,Aks(:,row_no,col_no),nvtx,nvtx);
            A;
        end
        b = b+sparse(nrow,1,bks(:,row_no),nvtx,1);
    end
    % get discrete Dirichlet boundary data
    b_nodes=find((xv==0)|(xv==1));
    size(b_nodes);
    int_nodes=1:nvtx;
    int_nodes(b_nodes)=[]; 
    b_int=b(int_nodes);
    wB=feval('g_eval',ns,xv(b_nodes),yv(b_nodes));
    % solve linear system for interior nodes;
    n_int = (ns+1)*(ns-1);
    A_int=A(int_nodes,int_nodes);
    A_ib =A(int_nodes,b_nodes);
    rhs=(b_int-A(int_nodes,b_nodes)*wB);
    u_int=A_int\rhs;
    A_int = reshape(A_int,n_int*n_int,1);
    A_int = full(A_int);
    uh=zeros(nvtx,1); uh(int_nodes)=u_int; uh(b_nodes)=wB;
    %size(uh)
    m = ns +1;
    uh = reshape(uh,m,m);
    uh = uh';
    uh = reshape(uh,nvtx,1);
end
%Dirichlet Data
function g=g_eval(ns,x,y)
    g=zeros(size(x));
    n_bound = (ns+1)*2;
    for i=1:2:n_bound
        g(i) = 10;
        g(i+1) =0;
    end    
end