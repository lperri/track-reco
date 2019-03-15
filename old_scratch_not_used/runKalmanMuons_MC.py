import FWCore.ParameterSet.Config as cms
process = cms.Process("L1KMTF4")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')   



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
"/store/mc/PhaseIIFall17D/DisplacedSUSY_SmuonToMuNeutralino_M-200_CTau-1000_TuneCUETP8M1_14TeV-pythia8/GEN-SIM-DIGI-RAW/L1TnoPU_93X_upgrade2023_realistic_v5-v1/80000/04B4FA4A-14A4-E811-991E-A0369FD0B17A.root" 
    ),
#    skipEvents=cms.untracked.uint32()
)








process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v5', '')

#Emulator Parameter
#process.load('L1Trigger.L1TMuonBarrel.fakeBmtfParams_cff')
process.load('L1Trigger.L1TMuonBarrel.simKBmtfStubs_cfi')
process.load('L1Trigger.L1TTwinMux.fakeTwinMuxParams_cff')
process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(
      cms.PSet(record = cms.string('L1TMuonBarrelParamsRcd'),
               data = cms.vstring('L1TMuonBarrelParams'))
                   ),
   verbose = cms.untracked.bool(True)
)
#process.fakeBmtfParams.fwVersion = cms.uint32(2)
#process.fakeBmtfParams.BX_max = cms.int32(2)
#process.fakeBmtfParams.BX_min = cms.int32(-2)

#enable stations here-each digit corresponds to a sector
#maskenable      = '000000000000'
#maskdisable     = '111111111111'
#process.fakeBmtfParams.mask_phtf_st1        = cms.vstring(maskdisable, maskenable, maskenable, maskenable, maskenable, maskenable, maskdisable)
#process.fakeBmtfParams.mask_phtf_st2        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
#process.fakeBmtfParams.mask_phtf_st3        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
#process.fakeBmtfParams.mask_phtf_st4        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)

#process.fakeBmtfParams.mask_ettf_st1        = cms.vstring(maskdisable, maskenable, maskenable, maskenable, maskenable, maskenable, maskdisable)
#process.fakeBmtfParams.mask_ettf_st2        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
#process.fakeBmtfParams.mask_ettf_st3        = cms.vstring(maskenable,  maskenable, maskenable, maskenable, maskenable, maskenable, maskenable)
####BMTF Emulator is loaded here
process.load('L1Trigger.L1TTwinMux.simTwinMuxDigis_cfi')
process.load('L1Trigger.L1TMuonBarrel.simBmtfDigis_cfi')
process.load('L1Trigger.L1TMuonBarrel.simKBmtfDigis_cfi')


process.simKBmtfStubs.srcTheta = cms.InputTag('simTwinMuxDigis')
process.simKBmtfStubs.srcPhi = cms.InputTag('simTwinMuxDigis')

process.simBmtfDigis.Debug = cms.untracked.int32(0)
process.simBmtfDigis.DTDigi_Source = cms.InputTag("simTwinMuxDigis")
process.simBmtfDigis.DTDigi_Theta_Source = cms.InputTag("simTwinMuxDigis")


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('displacedMuon140.root'),
    outputCommands = cms.untracked.vstring("drop *_*_*_*",
                                 "keep *_*_*_L1KMTF*",
                                 "keep *_muons_*_*",
                                 "keep *_simTwinMuxDigis_*_*",
                                 "keep *_simDtTrigger*_*_*",
                                 "keep *_dtTrigger*_*_*",
                                 "keep *_simBmtfDigis_*_*",
                                 "keep *_g4SimHits_MuonDTHits_*",
                                 "keep *_genParticles_*_*")                                       
)

  
#process.p = cms.Path(process.l1KalmanMuons)
process.p = cms.Path(process.simBmtfDigis*process.simKBmtfStubs*process.simKBmtfDigis)
process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.p,process.e)

