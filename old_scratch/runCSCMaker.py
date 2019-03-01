import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/scratch3/MTF/data/181025/singleMu140_PhaseII.root'
    )
)

process.cscTriggerHits = cms.EDProducer('CSCMaker'
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('cscHits.root')
)

  
process.p = cms.Path(process.cscTriggerHits)

process.e = cms.EndPath(process.out)
