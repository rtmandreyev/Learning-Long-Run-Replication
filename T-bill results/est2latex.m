% (c) Ida Johnsson, 2016 ida.b.johnsson@gmail.com

function f=est2latex(est,est_ref,rownames,colnames,caption,filename)

% Outputs latex table for m regression models with k coefficients
% estimated in each model.

% To adjust significance levels and *, **, ***, modify lines 69, 72 & 75
% To delete or modify notes under table, modify lines 128-131

% Make sure to either include the booktabs package in your latex file or
% replace \\toprule and \\midrule with \\
% Also, you need to include the package threeparttable

% Denote the j'th coefficient in model i by coef_ji,
% denote the stadard error of coef_ji by se_ji

%%      INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% est - (2*k) x m table of estimated coefficients and standard errors 
% formatted in the following way
%
% for i=1:m 
% for j=1:2:(2k-1)
% est(j,i)= coef_ji
% est(j+1,i)=se_ji
% end
% end   
%
% i.e.
%
% [coef_11 coef_12 ... coef_1m;
%  se_11   se_12   ... se_1m;
%       ....
%  coef_k1 coef_k2 ... coef_km;
%  se_k1   se_k2   ... se_km]

% rownames - 1 x k cell array with rownames, for example
% rownames = {'name_1', ..... , 'name_k'}

% colnames - 1 x m cell arra with column names, for example
% colnames = {'name_1', ... , 'name_m'}

% Note that latex symbols can be included, for example
% rownames = {'$\\beta_1$', .... , '$\\beta_k$'}

% caption - string containing the caption for the table, for example
% caption = 'Estimates of $\\mathbf{\\beta}$ for models $1$, $2$,
%            $\\ldots$, $m$'

% filename - name of latex file for storing the table, for example
% filename = 'table1.tex'


%% ADD SIGNIFICANCE STARS

[k2,m]=size(est);
k=k2/2;

% replace s.e. with t-stats

if isempty(est_ref)
    est_ref = zeros(k,m);
end

est_tst=est;

for i=2:2:2*k
est_tst(i,:)=(est(i-1,:)-est_ref(i/2,:))./est(i,:);   
     
end

stars=cell(2*k,m);

for i=2:2:(2*k)
    for j=1:m
    
        if abs(est_tst(i,j))>1.65 && abs(est_tst(i,j))<=1.96
            stars(i-1,j)={'*'};
            
        elseif abs(est_tst(i,j))>1.96 && abs(est_tst(i,j))<=2.58
            stars(i-1,j)={'**'};
            
        elseif abs(est_tst(i,j))>2.58
            stars(i-1,j)={'***'};
        
        end
        
    end
        
end 

%% CREATE LATEX TABLE

[rows,m]=size(est);
k=rows/2;


FID = fopen(filename, 'w');
fprintf(FID, '\\begin{table}[!h]\\caption{%s}\\centering \n',caption);
fprintf(FID, '\\begin{threeparttable} \n');

fprintf(FID, '\\begin{tabular}{|l|*{%.0f}{c}|}\\toprule \n',m);

% column names
for j=1:m
    fprintf(FID,'& %s',colnames{j});
end

fprintf(FID,'\\\\ \\midrule \n');

% 

for i=1:(2*k) % rows
    
     if mod(i,2) == 1 % odd row - estimated coefficients
fprintf(FID, '%s& %.2f%s ',['\multirow{2}{*}{' rownames{(i-1)/2+1} '}'],est(i,1), stars{i,1});
     else % even row - s.e.
fprintf(FID, '& (%.2f) ',est(i,1));
     end
    
    
    for j=2:(m-1) % columns
        if mod(i,2) == 1 % odd row - estimated coefficients
fprintf(FID, '& %.2f%s ',est(i,j), stars{i,j});
        else % even row - s.e.
fprintf(FID, '& (%.2f) ',est(i,j));
        end
    end
    
    if mod(i,2) == 1 % last column
fprintf(FID, '& %.2f%s \\\\ \n',est(i,m), stars{i,m});
    else
        
fprintf(FID, '& (%.2f) \\\\ \n',est(i,m));
    end
    
end


fprintf(FID,'\\bottomrule \\end{tabular} \n');

%tablenotes
fprintf(FID,'\\begin{tablenotes}[para, flushleft] \n');
fprintf(FID,'\\item Standard errors are in parentheses, *, ** and *** denote statistical significance at 10 percent, 5 percent and 1 percent levels, respectively. \\\\ \n');
fprintf(FID,'\\item Note 2 \n');
fprintf(FID,'\\end{tablenotes} \n');
fprintf(FID,'\\end{threeparttable} \n');
fprintf(FID,'\\end{table} \n');
fclose(FID);
    



end