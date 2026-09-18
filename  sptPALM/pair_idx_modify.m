% make the reference chart to map the cell idx
% #17 
pairchart12 = pairchart;
pairchart12{2,6} = [1 2 4 3];
pairchart12{2,7} = [1 4 2 5 7 6 1000];
pairchart12{4,8} = [2 3 1 4]; 
save('pairchart',"pairchart12");

%%
% 42
pairchart12 = pairchart;
pairchart12{2,1} = [2 1 3 4];
pairchart12{1,2} = [1 2 3 4];
pairchart12{2,2} = [100 1 2 3];
save('pairchart',"pairchart12");