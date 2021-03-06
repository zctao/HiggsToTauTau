import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.Framework.test.cmsExceptionsFatal_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")

process.source = cms.Source("EmptySource",
     firstLuminosityBlock = cms.untracked.uint32(1),
#     numberEventsInLuminosityBlock = cms.untracked.uint32(200),
     firstEvent = cms.untracked.uint32(1),
     firstRun = cms.untracked.uint32(1),
     numberEventsInRun = cms.untracked.uint32(10000)

)

from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(8000.),
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
            TauolaPolar,
                TauolaDefaultInputCards
        ),
        parameterSets = cms.vstring('Tauola')
        ), 
        
    PythiaParameters = cms.PSet(
        py8ProcessSettings = cms.vstring( 'StringZ:usePetersonB = on', # these 2 together ==
                                          'StringZ:usePetersonC = on', # mstj(11)=3
                                          'HiggsSM:all = on' # 
                                          # what about choice of structure function ??? (mstp(51)=7)
                                          ),
        py8HDecaySettings = cms.vstring(  '25:onMode = off', # turn OFF all Higgs decays
                                          '25:onIfAny = 15'  # turn ON H->tautau
                                          ),
        py8TauDecaySettings = cms.vstring('15:onMode = off', # turn OFF all tau decays
                                          ),
        parameterSets = cms.vstring(
                                          'py8ProcessSettings',
                                          'py8HDecaySettings')
        )
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger",
#    cout = cms.untracked.PSet(
#        default = cms.untracked.PSet(
#            limit = cms.untracked.int32(10000)
#        )
#    ),
#    destinations = cms.untracked.vstring('cout')
#)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    generator = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('HepJamesRandom')
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.GEN = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('mcSamples/8TeV/all/Py8_H2tautau_tauola.root')
)

process.printGenParticles = cms.EDAnalyzer("ParticleListDrawer",
                                           src = cms.InputTag("genParticles"),
                                           maxEventsToPrint = cms.untracked.int32(3) )

process.p = cms.Path(process.generator*process.genParticles*process.printGenParticles)
process.outpath = cms.EndPath(process.GEN)

process.schedule = cms.Schedule(process.p, process.outpath)
