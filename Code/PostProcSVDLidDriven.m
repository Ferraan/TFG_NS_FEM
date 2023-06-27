%%Compute all subspace calculations for all possible permutations
%As long as subspace(A)>=subspace(B)
disp("Calculating angles between subspaces LidDriven")
disp("********************************************")
disp("********************************************")
files=dir("GIDPOST/LidDriven75NSS*");
combinations=nchoosek(files,2);
rows=1;
columns=size(combinations,1);
for i=1:size(combinations,1)
    load(strcat("GIDPOST\",combinations(i,1).name,'/',combinations(i,1).name,".mat"));
    U1=SNAP_cluster.u.U;
    load(strcat("GIDPOST\",combinations(i,2).name,'/',combinations(i,2).name,".mat"));
    U2=SNAP_cluster.u.U;
    if(max(size(U1,2),size(U2,2))>rows)
        rows=max(size(U1,2),size(U2,2))+2;
    end
    if(size(U1,2)>=size(U2,2)) %p>=q
        [Q1,R1]=qr(U1,"econ");
        [Q2,R2]=qr(U2,"econ");
        C=Q1.'*Q2;
        [Y,SIGMA,Z]=svd(C);
        cosangles=diag(SIGMA);
        disp(strcat("Angles between ",combinations(i,1).name," and ",combinations(i,2).name))
        disp(strcat('Dimensions of first subspace: ',num2str(size(U1,2)),' Dimensions of second subspace: ',num2str(size(U2,2))))
        disp(cosangles);
        %%%FIND Re and steps
        %Re
        pat="_U_";
        str=combinations(i,1).name;
        pat2="_mu_";
        TF2=strfind(str,pat2);
        TF=strfind(str,pat);
        mu1=str2double(str(TF2+4:TF-1));
        Re1=1/mu1;
        %Steps
        patSteps="Steps_";
        patStepsEND="_tolU";
        TF=strfind(str,patSteps);
        TF2=strfind(str,patStepsEND);
        steps1=str2double(str(TF+6:TF2-1));
        
        %%%%%
        str=combinations(i,2).name;
        TF2=strfind(str,pat2);
        TF=strfind(str,pat);
        mu2=str2double(str(TF2+4:TF-1));
        Re2=1/mu2;
        patSteps="Steps_";
        patStepsEND="_tolU";
        TF=strfind(str,patSteps);
        TF2=strfind(str,patStepsEND);
        steps2=str2double(str(TF+6:TF2-1));
        disp(strcat(num2str(Re1)," ",num2str(steps1)," ",num2str(Re2)," ",num2str(steps2)));
    else   
       
        U3=U1;
        U1=U2;
        U2=U3;
        [Q1,R1]=qr(U1,"econ");
        [Q2,R2]=qr(U2,"econ");
        C=Q1.'*Q2;
        [Y,SIGMA,Z]=svd(C);
        cosangles=diag(SIGMA);
        disp(strcat("Angles between ",combinations(i,2).name," and ",combinations(i,1).name))
        disp(strcat('Dimensions of first subspace: ',num2str(size(U1,2)),' Dimensions of second subspace: ',num2str(size(U2,2))))
        disp(cosangles);
        %%%FIND Re and steps
        %Re
        pat="_U_";
        str=combinations(i,2).name;
        pat2="_mu_";
        TF2=strfind(str,pat2);
        TF=strfind(str,pat);
        mu1=str2double(str(TF2+4:TF-1));
        Re1=1/mu1;
        %Steps
        patSteps="Steps_";
        patStepsEND="_tolU";
        TF=strfind(str,patSteps);
        TF2=strfind(str,patStepsEND);
        steps1=str2double(str(TF+6:TF2-1));
        
        %%%%%
        str=combinations(i,1).name;
        TF2=strfind(str,pat2);
        TF=strfind(str,pat);
        mu2=str2double(str(TF2+4:TF-1));
        Re2=1/mu2;
        patSteps="Steps_";
        patStepsEND="_tolU";
        TF=strfind(str,patSteps);
        TF2=strfind(str,patStepsEND);
        steps2=str2double(str(TF+6:TF2-1));
        disp(strcat(num2str(Re1)," ",num2str(steps1)," ",num2str(Re2)," ",num2str(steps2)));

    end
    Reynoldssup(i)=Re1;
    Reynoldsinf(i)=Re2;
    Stepssup(i)=steps1;
    Stepsinf(i)=steps2;
    angles(1:size(cosangles,1),i)=cosangles;
end
%%
%%print matrix for latex
disp("*****************************")
disp("MATRIX LATEX:")
disp("*****************************")
ntables=4;

begincol=1;
columns=floor(columns/ntables)+1;
for k=1:4
    for i=1:rows
        for j=begincol:columns
            if(i==1)
                fprintf('Re=%d \\(\\nsteps\\)=%d', Reynoldssup(j),Stepssup(j));
                
            elseif(i==2)
                fprintf('Re=%d \\(\\nsteps\\)=%d', Reynoldsinf(j),Stepsinf(j));
            else
                if(angles(i-2,j)~=0)
                fprintf("%1.4f",angles(i-2,j))
                else
                    fprintf("-")
                end
            end
            if(j~=columns)
                fprintf(" & ")
            end
        end
        disp("\\ \hline")
    end
    begincol=columns+1;
    columns=columns+floor(15/ntables);
    disp("NEXT****************")
end
%%