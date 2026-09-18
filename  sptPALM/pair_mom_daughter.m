clear;close all;
 % cell type 1 means cell type 1/2, which are files with 3/4 suffix; 
% 3 means cell type 3/4, which are files with 5/6 suffix; 

% find the linked cell labels 
dataall = importdata('42_20201211_merge_Data2Plot_thre0_RG_Z_01-7.mat');
pairchart = importdata('pairchart.mat');
cell1reform = reformatrix(dataall(1,:));
cell2reform = reformatrix(dataall(2,:));
cell3reform = reformatrix(dataall(3,:));
cell4reform = reformatrix(dataall(4,:));
% pair the matrix by using pairchart of the dataset 
[momnew1,daunew1] = reformthroghtchart(cell1reform,cell2reform,pairchart,1);
[momnew3,daunew3] = reformthroghtchart(cell3reform,cell4reform,pairchart,3);
%%
[daughtersizecor2,daughtersizecor2fra] = correctmothercellsize(momnew1,daunew1,100,10,0.45);
[daughtersizecor2_2,daughtersizecor2fra_2] = correctmothercellsize(momnew3,daunew3,100,10,0.45);
figure
plot(daughtersizecor2,daughtersizecor2fra,'r.','MarkerSize',17);
axis([min(daughtersizecor2)-5 max(daughtersizecor2)+5 0 1]);
hold on 
plot(daughtersizecor2_2,daughtersizecor2fra_2,'g.','MarkerSize',17);
%%
% find the mother cell size below certain threshold 
function [daughtersizecor2,daughtersizecor2fra] = correctmothercellsize (momnew1,daunew1,cellvolthre,tracknumthre,scalefactor)
smallmothercell = find(momnew1(:,5)<cellvolthre); % cell volume threshold; 
daughtersmallmom = daunew1(smallmothercell,:);
smallmothercell2 = momnew1(smallmothercell,:);

daughtersmallmom2 = find(daughtersmallmom(:,7)>tracknumthre); % track number threshold 
daughtersmallmom3 = daughtersmallmom(daughtersmallmom2,:);
%%
% correct the daughter cell size with mother cell size, with the factor
% 0.45; 
daughtersizecoor = daughtersmallmom(:,5) - smallmothercell2(:,5).*scalefactor; 
daughtersizecor2 = daughtersizecoor(daughtersmallmom2,:);
daughtersizecor2fra = daughtersmallmom3(:,6);

end

%%
function [momnew,daunew] = reformthroghtchart(cell1reform,cell2reform,pairchart,celltype)
count = 1;
for momcelliter = 1:size(cell1reform,1)
    % find the cell index 
    momfieldidx = cell1reform(momcelliter,1);
    momcellidx = cell1reform(momcelliter,2);
    if celltype == 1
    daucellidx = pairchart{2,momfieldidx};
    else 
        daucellidx = pairchart{4,momfieldidx};
    end
    daucellidx2 = daucellidx(momcellidx);
    daughternewidx = find(cell2reform(:,1) == momfieldidx & cell2reform(:,2) == daucellidx2);
    if isempty(daughternewidx)
        continue
    else 
        momnew(count,:) = cell1reform(momcelliter,:);
        daunew (count,:) = cell2reform(daughternewidx,:);
        count = count +1;
    end
end
% convert the cell vol to fL 
momnew(:,5) = momnew(:,5).*0.16^3;
daunew(:,5) = daunew(:,5).*0.16^3;
momnew(:,4) = momnew(:,4).*0.16^2;
daunew(:,4) = daunew(:,4).*0.16^2;
end

%%
% function to make the new matrix 
function datareform = reformatrix(datacelltype1)
BFidx = datacelltype1{:,6};
BFfra = datacelltype1{:,2};
BFnum = datacelltype1{:,4};
datareform = [BFidx,BFfra,BFnum];
end
