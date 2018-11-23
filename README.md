# ISR2016MuonUnfold

make lib

root
.L ISR2016MuonUnfold.cc
SaveBinning()
SaveHist()
Unfold()
Compare("hdataUnfold","hmcGen_GenAxis",0,80,100)

