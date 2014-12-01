import FWCore.ParameterSet.Config as cms

process = cms.Process("ALL")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# to run over the generated sample
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/z/ztao/public/CMSSW_6_2_0_SLHC12_patch1/src/Configuration/Generator/python/mcSamples/13TeV/gg2H/Py8_H2tautau_tauola.root')
)



# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')




# ---------------------------------------------------------------------------
#
#process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')


process.GeneratorNtupleMaker = cms.EDAnalyzer("GeneratorNtupleMaker",
                                              MyProcess = cms.int32(25)
)

process.genH = cms.Path(process.GeneratorNtupleMaker)


# ---------------------------------------------------------------------------



process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('RootFiles/13TeV/H2tautauNtuple_gg2H.root'),
)
