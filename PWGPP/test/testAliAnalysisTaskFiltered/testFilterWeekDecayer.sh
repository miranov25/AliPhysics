# source $AliPhysics_SRC/PWGPP/test/testAliAnalysisTaskFiltered/testFilterWeekDecayer.sh

alias helpCat=cat
[[ -x "$(command -v pygmentize)" ]] && alias helpCat="pygmentize -O style=borland,linenos=1 -l bash"

init(){
    [[ -z "${ALILOG_HOST}" ]] && source ${ALICE_ROOT}/libexec/alilog4bash.sh
    [[ -z "${PWGPP_runMap}" ]] && source ${ALICE_ROOT}/libexec/utilities.sh

  cat <<HELP_USAGE | helpCat
    # runParallel
    init
    makeESDList
    splitESDlist
    runParallel
    #
    makeNDMapsInvMassParallel
HELP_USAGE

}

makeESDList(){
  # dummy function to create a list
  # below test MC perfromance productiion PbPb hijing + flat 1pt jets +  350 chunks x5 events
  for esd in $(find /lustre/alice/DETdata/alice/sim/2020/LHC20f9a/ -iname "root_archive.zip" | grep -v AODFilterTrees) ; do echo $esd#AliESDs.root; done > esd.list
}

splitESDlist(){
    cat <<HELP_USAGE | helpCat
    splitESDlist
    Parameters:
        1 - n chuks to split
    Algorithm:
HELP_USAGE
   nChunks=$1
   split --suffix-length=3 -d -l ${nChunks} esd.list dir
   for a in $(ls -d dir*); do
     echo $a
     mkdir -p ${a}Split
     mv ${a}  ${a}Split/esd.list
   done
   ls -d  dir*Split* > dir.list
}

runParallel(){
  cat <<HELP_USAGE | helpCat
    runParallel
    Parameters:
        1 - nChunks
        2 - nEvents
        3 - nJobs
        4 -
    Algorithm:
HELP_USAGE
    [[ $# -ne 3 ]] &&return
    export nEvents=$2
    export nChunks=$1
    export nJobs=$3
    cat <<EOF >  runParallel.sh
#!/bin/bash
    alilog_info "runParallel Begin" $(date)
    touch Begin.info
    root.exe -n -b -l <<\EOF 2>&1 | tee runParalel.log
    .L $AliPhysics_SRC/PWGPP/test/testAliAnalysisTaskFiltered/AliAnalysisTaskFilteredTest.C
    // TODO OCDB as shell parameter
    TStopwatcher timer;
    //AliAnalysisTaskFilteredTest("esd.list",0,1,1,1,"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2018/OCDB/",${nChunks},0,${nEvents},0,1,1);
    //AliAnalysisTaskFilteredTest("esd.list",0,1,1,1,"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2018/OCDB/",${nChunks},0,${nEvents},0,1,0);
    AliAnalysisTaskFilteredTest("esd.list",0,1000,100,10000,"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2018/OCDB/",${nChunks},0,${nEvents},0,1,1);
     timer.Print();
    .q
    alilog_info "runParallel END" $(date)
    touch end.info
EOF
   chmod a+x runParallel.sh
   export mDir=$(pwd)
   head -n ${nChunks} dir.list   | parallel --memfree 4G -j${nJobs} " cd {};  ${mDir}/runParallel.sh  >  tee runParallel.log"
    for a in $(find  dir*/ -size +1k  -iname  Filter*.root | xargs dirname); do echo $a/AliAnalysisTaskWeakDecayVertexer.root;done  > debug.list
   rm  AliAnalysisTaskWeakDecayVertexer.root
   alihadd   -k AliAnalysisTaskWeakDecayVertexer.root dir*/AliAnalysisTaskWeakDecayVertexer.root
   ali_info runParallelxxx
}

testMerge(){
   # this is dummy test file
    rm AliAnalysisTaskWeakDecayVertexer.root

    #
    for a in $(find  dir*/ -size +1k  -iname  Filter*.root | xargs dirname); do echo $a/AliAnalysisTaskWeakDecayVertexer.root;done  > debug.list
    alihadd -k AliAnalysisTaskWeakDecayVertexer.root @debug.list

}

makeNDMapsInvMassParallel(){
 cat <<HELP_USAGE | helpCat
    makeNDMaps - make maps for n=8 histograms - if more histograms added number has to be modified
    it takes ~ 10 minutes to extract maps
    Parameters:
      $1 - nJobs
    Algorithm:
HELP_USAGE
    [[ $# -ne 1 ]] &&return
    cat <<EOF >  makendPipelineInvMassl.sh
#!/bin/bash
    export hisIndex=\$1
    root.exe -n -b -l <<\EOF 2>&1 | tee ndPipelineInvMass.log
    .L $NOTES/JIRA/ATO-544/analyzeHybrid.C
    int iMap= TString(gSystem->Getenv("hisIndex")).Atoi();
    makeMaps(iMap);
    .q
EOF
   chmod a+x  makendPipelineInvMassl.sh
   seq 0 8 |  parallel --memfree 4G -j"${nJobs}" " ./makendPipelineInvMassl.sh {} | tee makeNDMapInvMass{}.log"
   rm -f mapInvariantMass.root
   alihadd -k mapInvariantMass.root mapInvariantMass_*.root
}



init