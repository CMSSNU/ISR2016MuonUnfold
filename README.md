# ISR2016MuonUnfold

make lib

root  
.L ISR2016MuonUnfold.cc  
SaveHist("hist.root")  
Unfold("hist.root","unfold.root",MODE,SCANMETHOD,f)  
//MODE -1:PrivateRegularization 0:NoRegularization 1:SizeRegularization 2:Derivative 3:Curvacure  
 
Compare("unfold.root","hdataUnfold","hmcGen_GenAxis",0,80,100)  

