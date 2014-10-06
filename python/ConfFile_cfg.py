import FWCore.ParameterSet.Config as cms

process = cms.Process("OwnParticles")

process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometryPostLS1XML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load('Configuration/StandardSequences/GeometryIdeal_cff')


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("MuonTools.StandardSequences.RecoStandAloneMuon_cff")

# have stareco use hlt digis
process.load("MuonTools.Configuration.HltMuonDigis_cff")
# have stareco use hlt segments (1)
process.load("MuonTools.Configuration.HltMuonSegments_cff")
# keep stareco from using rpcs
process.load("MuonTools.Configuration.StandAloneNoRpc_cff")

process.GlobalTag.globaltag = 'GR_P_V32::All'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/d/depoyraz/RPC/sample209151.root'
    )
)

process.dTandCSCSegmentsinTracks = cms.EDProducer("DTandCSCSegmentsinTracks",
                                                  cscSegments = cms.InputTag("hltCscSegments"),
                                                  dt4DSegments = cms.InputTag("hltDt4DSegments"),
                                                  tracks = cms.InputTag("standAloneMuons","UpdatedAtVtx")
                                                  )

process.rpcPointProducer = cms.EDProducer('RPCPointProducer',
  incldt = cms.untracked.bool(True),
  inclcsc = cms.untracked.bool(True),
  incltrack =  cms.untracked.bool(False),

  debug = cms.untracked.bool(False),

  rangestrips = cms.untracked.double(4.),
  rangestripsRB4 = cms.untracked.double(4.),
  MinCosAng = cms.untracked.double(0.85),
  MaxD = cms.untracked.double(80.0),
  MaxDrb4 = cms.untracked.double(150.0),
  ExtrapolatedRegion = cms.untracked.double(0.6), #in stripl/2 in Y and stripw*nstrips/2 in X
  cscSegments = cms.InputTag('dTandCSCSegmentsinTracks','SelectedCscSegments','OwnParticles'),
  dt4DSegments = cms.InputTag('dTandCSCSegmentsinTracks','SelectedDtSegments','OwnParticles'),
  tracks = cms.InputTag("standAloneMuons"),
  TrackTransformer = cms.PSet(
      DoPredictionsOnly = cms.bool(False),
      Fitter = cms.string('KFFitterForRefitInsideOut'),
      TrackerRecHitBuilder = cms.string('WithTrackAngle'),
      Smoother = cms.string('KFSmootherForRefitInsideOut'),
      MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
      RefitDirection = cms.string('alongMomentum'),
      RefitRPCHits = cms.bool(False),
      Propagator = cms.string('SmartPropagatorAnyRKOpposite')
  )
)


process.panalyzer = cms.EDAnalyzer('RPCPointAnalyzer',
 							   RootFileName = cms.untracked.string("PointAnalyzerHistograms.root"),
							   #rpcDTPoints = cms.InputTag("rpcPointProducer","RPCDTExtrapolatedPoints"),
							   #rpcCSCPoints = cms.InputTag("rpcPointProducer","RPCCSCExtrapolatedPoints"),
							   rpcDTPoints = cms.InputTag("rpcPointProducer","RPCDTExtrapolatedPoints"),
							   #rpcDTPoints = cms.InputTag("hltRPCDTExtrapolatedPoints"),
  							   rpcRecHits = cms.InputTag("hltRpcRecHits"),


)


process.normfilter = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring("AlCa_RPCMuonNormalisation*"),
    eventSetupPathsKey = cms.string(''),
    andOr = cms.bool(True),
    throw = cms.bool(True)
)

process.p = cms.Path(process.normfilter*process.muonstandalonereco*process.dTandCSCSegmentsinTracks*process.rpcPointProducer*process.panalyzer)
