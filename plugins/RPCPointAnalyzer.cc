// -*- C++ -*-
//
// Package:    RPCAnalyzer/RPCPointAnalyzer
// Class:      RPCPointAnalyzer
// 
/**\class RPCPointAnalyzer RPCPointAnalyzer.cc RPCAnalyzer/RPCPointAnalyzer/plugins/RPCPointAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Deniz Poyraz
//         Created:  Fri, 18 Jul 2014 09:33:03 GMT
//
//


// system include files
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


// root include files
#include <TRandom.h>
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///Data Format
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>

///Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>
///Log messages
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoLocalMuon/RPCRecHit/interface/RPCPointProducer.h"

//
// class declaration
//

class RPCPointAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RPCPointAnalyzer(const edm::ParameterSet&);
      ~RPCPointAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::ESHandle<RPCGeometry> rpcGeom;
      //edm::ESHandle<DTGeometry> dtGeo;
      //edm::ESHandle<CSCGeometry> cscGeo;
      edm::InputTag rpcRecHitLabel;
      edm::InputTag rpcDTPointsLabel;
      edm::InputTag rpcCSCPointsLabel;
 	 // std::string digiLabel, rootFileName;
 	  std::string  rootFileName;
  	  TFile * outputfile;
  	  
  	  TH1F * MB1_R0_Sec1_Lay1_Local_RecHits,	*MB1_R0_Sec1_Lay1_Local_RecHits_Error,  *MB1_R0_Sec1_Lay1_Local_DTPoints,	*MB1_R0_Sec1_Lay1_Local_DTPoints_Error,  * MB1_R0_Sec1_Lay1_Local_Difference;
  	  TH1F * MB1_R0_Sec1_Lay2_Local_RecHits,	*MB1_R0_Sec1_Lay2_Local_RecHits_Error,  *MB1_R0_Sec1_Lay2_Local_DTPoints,	*MB1_R0_Sec1_Lay2_Local_DTPoints_Error,  * MB1_R0_Sec1_Lay2_Local_Difference;
  	   
  	   
  	  TH1F * MB1_R0_Sec2_Lay1_Local_RecHits,	*MB1_R0_Sec2_Lay1_Local_RecHits_Error,  *MB1_R0_Sec2_Lay1_Local_DTPoints,	*MB1_R0_Sec2_Lay1_Local_DTPoints_Error,  * MB1_R0_Sec2_Lay1_Local_Difference;
  	  TH1F * MB1_R0_Sec2_Lay2_Local_RecHits,	*MB1_R0_Sec2_Lay2_Local_RecHits_Error,  *MB1_R0_Sec2_Lay2_Local_DTPoints,	*MB1_R0_Sec2_Lay2_Local_DTPoints_Error,  * MB1_R0_Sec2_Lay2_Local_Difference;
      
      int n_x  = 100;  double n1_x  = -150,  n2_x  = 150;
      int n_dx  = 50;  double n1_dx  = 0,  n2_dx  = 5;
      int n_ddx  = 100;  double n1_ddx  = 1,  n2_ddx  = 2;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RPCPointAnalyzer::RPCPointAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   
   //****  also need to be defined in the py code
   
   rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
   rpcRecHitLabel     =  iConfig.getParameter<edm::InputTag>("rpcRecHits");
   rpcDTPointsLabel  = iConfig.getParameter<edm::InputTag>("rpcDTPoints");
   rpcCSCPointsLabel  = iConfig.getParameter<edm::InputTag>("rpcCSCPoints");
	
	outputfile = new TFile(rootFileName.c_str(), "RECREATE" );
   
    MB1_R0_Sec1_Lay1_Local_RecHits = new TH1F("MB1_R0_Sec1_Lay1_RecHits","MB1_R0_Sec1_Lay1_RecHits",n_x,n1_x,n2_x);
    MB1_R0_Sec1_Lay1_Local_RecHits_Error = new TH1F("MB1_R0_Sec1_Lay1_RecHits_Error","MB1_R0_Sec1_Lay1_RecHits_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec1_Lay1_Local_DTPoints = new TH1F("MB1_R0_Sec1_Lay1_DTPoints","MB1_R0_Sec1_Lay1_DTPoints",n_x,n1_x,n2_x);
    MB1_R0_Sec1_Lay1_Local_DTPoints_Error = new TH1F("MB1_R0_Sec1_Lay1_DTPoints_Error","MB1_R0_Sec1_Lay1_DTPoints_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec1_Lay1_Local_Difference = new TH1F("MB1_R0_Sec1_Lay1_Difference","MB1_R0_Sec1_Lay1_Difference",n_dx,n1_dx,n2_dx);
    
    MB1_R0_Sec1_Lay2_Local_RecHits = new TH1F("MB1_R0_Sec1_Lay2_Local_RecHits","MB1_R0_Sec1_Lay2_Local_RecHits",n_x,n1_x,n2_x);
    MB1_R0_Sec1_Lay2_Local_RecHits_Error = new TH1F("MB1_R0_Sec1_Lay2_Local_RecHits_Error","MB1_R0_Sec1_Lay2_Local_RecHits_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec1_Lay2_Local_DTPoints = new TH1F("MB1_R0_Sec1_Lay2_Local_DTPoints","MB1_R0_Sec1_Lay2_Local_DTPoints",n_x,n1_x,n2_x);
    MB1_R0_Sec1_Lay2_Local_DTPoints_Error = new TH1F("MB1_R0_Sec1_Lay2_Local_DTPoints_Error","MB1_R0_Sec1_Lay2_Local_DTPoints_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec1_Lay2_Local_Difference = new TH1F("MB1_R0_Sec1_Lay2_Local_Difference","MB1_R0_Sec1_Lay2_Local_Difference",n_dx,n1_dx,n2_dx);
    
    
    MB1_R0_Sec2_Lay1_Local_RecHits = new TH1F("MB1_R0_Sec2_Lay1_RecHits","MB1_R0_Sec2_Lay1_RecHits",n_x,n1_x,n2_x);
    MB1_R0_Sec2_Lay1_Local_RecHits_Error = new TH1F("MB1_R0_Sec2_Lay1_RecHits_Error","MB1_R0_Sec2_Lay1_RecHits_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec2_Lay1_Local_DTPoints = new TH1F("MB1_R0_Sec2_Lay1_DTPoints","MB1_R0_Sec2_Lay1_DTPoints",n_x,n1_x,n2_x);
    MB1_R0_Sec2_Lay1_Local_DTPoints_Error = new TH1F("MB1_R0_Sec2_Lay1_DTPoints_Error","MB1_R0_Sec2_Lay1_DTPoints_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec2_Lay1_Local_Difference = new TH1F("MB1_R0_Sec2_Lay1_Difference","MB1_R0_Sec2_Lay1_Difference",n_dx,n1_dx,n2_dx);
    
    MB1_R0_Sec2_Lay2_Local_RecHits = new TH1F("MB1_R0_Sec2_Lay2_Local_RecHits","MB1_R0_Sec2_Lay2_Local_RecHits",n_x,n1_x,n2_x);
    MB1_R0_Sec2_Lay2_Local_RecHits_Error = new TH1F("MB1_R0_Sec2_Lay2_Local_RecHits_Error","MB1_R0_Sec2_Lay2_Local_RecHits_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec2_Lay2_Local_DTPoints = new TH1F("MB1_R0_Sec2_Lay2_Local_DTPoints","MB1_R0_Sec2_Lay2_Local_DTPoints",n_x,n1_x,n2_x);
    MB1_R0_Sec2_Lay2_Local_DTPoints_Error = new TH1F("MB1_R0_Sec2_Lay2_Local_DTPoints_Error","MB1_R0_Sec2_Lay2_Local_DTPoints_Error",n_dx,n1_dx,n2_dx);
    MB1_R0_Sec2_Lay2_Local_Difference = new TH1F("MB1_R0_Sec2_Lay2_Local_Difference","MB1_R0_Sec2_Lay2_Local_Difference",n_dx,n1_dx,n2_dx);
    
    

}


RPCPointAnalyzer::~RPCPointAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
   MB1_R0_Sec1_Lay1_Local_RecHits->Write();
   MB1_R0_Sec1_Lay1_Local_RecHits_Error->Write();
   MB1_R0_Sec1_Lay1_Local_DTPoints->Write();
   MB1_R0_Sec1_Lay1_Local_DTPoints_Error->Write();
   MB1_R0_Sec1_Lay1_Local_Difference->Write();
   
   MB1_R0_Sec1_Lay2_Local_RecHits->Write();
   MB1_R0_Sec1_Lay2_Local_RecHits_Error->Write();
   MB1_R0_Sec1_Lay2_Local_DTPoints->Write();
   MB1_R0_Sec1_Lay2_Local_DTPoints_Error->Write();
   MB1_R0_Sec1_Lay2_Local_Difference->Write();
   
   MB1_R0_Sec2_Lay1_Local_RecHits->Write();
   MB1_R0_Sec2_Lay1_Local_RecHits_Error->Write();
   MB1_R0_Sec2_Lay1_Local_DTPoints->Write();
   MB1_R0_Sec2_Lay1_Local_DTPoints_Error->Write();
   MB1_R0_Sec2_Lay1_Local_Difference->Write();
   
   MB1_R0_Sec2_Lay2_Local_RecHits->Write();
   MB1_R0_Sec2_Lay2_Local_RecHits_Error->Write();
   MB1_R0_Sec2_Lay2_Local_DTPoints->Write();
   MB1_R0_Sec2_Lay2_Local_DTPoints_Error->Write();
   MB1_R0_Sec2_Lay2_Local_Difference->Write();
   
   
}



//
// member functions
//

// ------------ method called for each event  ------------
void
RPCPointAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

	edm::Handle<RPCRecHitCollection> rpcRecHits; 
    iEvent.getByLabel(rpcRecHitLabel,rpcRecHits);
   
    edm::Handle<RPCRecHitCollection> DTPoints;
    iEvent.getByLabel(rpcDTPointsLabel,DTPoints);   
    
    edm::Handle<RPCRecHitCollection> CSCPoints;
    iEvent.getByLabel(rpcCSCPointsLabel,CSCPoints); 
    
    
  //********************	RPC RecHit Points	********************	    	 
	RPCRecHitCollection::const_iterator recHit;
	for(recHit = rpcRecHits->begin(); recHit != rpcRecHits->end(); recHit++){ //Filling DQMOccupancies
	 	 
		//int cls = recHit->clusterSize();
		//int firststrip = recHit->firstClusterStrip();
		//int bx = recHit->BunchX();
		RPCDetId rpcId = recHit->rpcId();
		LocalPoint rpc_position=recHit->localPosition();
		LocalError rpc_position_error=recHit->localPositionError();
		
		int region = rpcId.region();
		if(region == 0){
		const RPCRoll* rollasociated = rpcGeom->roll(rpcId);
		const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
		float stripw = top_->pitch();
		//float stripl = top_->stripLength();
		//float rpc_error_x = rpc_position_error.xx();
		float error_cal = stripw * stripw  /sqrt(12);
		const float stripp = rollasociated->strip(LocalPoint(rpc_position.x(),rpc_position.y(),0.));
		LocalError RollError = rollasociated->localError(stripp);
		

		std::cout<<" RPC Id: "<<rpcId<<std::endl;
		std::cout<<" RPC Position Error: "<<rpc_position_error<<std::endl;
		std::cout<<" RollError: "<<RollError<<" Error Calculated: "<<error_cal<<std::endl;
		}
			
	}
    
		    	 
//********************	DT Points	********************	    	 
//********************	Barrel		********************	    	 
		   	
	if(DTPoints.isValid() && DTPoints->begin()!=DTPoints->end()){ //No Empty Predictions
	RPCRecHitCollection::const_iterator rpcDTPoint;

		for(rpcDTPoint = DTPoints->begin(); rpcDTPoint != DTPoints->end(); rpcDTPoint++){	//DT Points
		LocalPoint PointExtrapolatedRPCFrame = rpcDTPoint->localPosition();
		LocalError PointExtrapolatedRPCFrameError = rpcDTPoint->localPositionError();
		double LocalDTPoint_x = PointExtrapolatedRPCFrame.x();
		double LocalDTPoint_x_error = PointExtrapolatedRPCFrameError.xx();
		RPCDetId  rpcId = rpcDTPoint->rpcId();

		//std::cout<<"RPC Id: "<<rpcId<<" PointExtrapolatedRPCFrame: "<<LocalDTPoint_x<<" Error: "<<LocalDTPoint_x_error<<std::endl;
		
		
		// this change
		
		typedef std::pair<RPCRecHitCollection::const_iterator, RPCRecHitCollection::const_iterator> rangeRecHits;
		rangeRecHits recHitCollection = rpcRecHits->get(rpcId);
		RPCRecHitCollection::const_iterator dtrecHit;

			for(dtrecHit = recHitCollection.first; dtrecHit != recHitCollection.second ; dtrecHit++) {
			//countRecHits++;
			LocalPoint recHitPos=dtrecHit->localPosition();
			LocalError recHitPosError=dtrecHit->localPositionError();
			double LocalRecHit_x = recHitPos.x();
			double LocalRecHit_x_error = recHitPosError.xx();
			double LocalRecHit_y = recHitPos.y();
			double LocalRecHit_y_error = recHitPosError.yy();
			double LocalRecHit_xy_error = recHitPosError.xy();
			
			double Local_diff_x = fabs(LocalRecHit_x - LocalDTPoint_x );
			int layer = rpcId.layer();
			int sector = rpcId.sector();
			
			const RPCRoll* rollasociated = rpcGeom->roll(rpcId);
			const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
			float stripw = top_->pitch();
			float error = stripw * stripw / sqrt(12);
			const float stripp = rollasociated->strip(LocalPoint(LocalDTPoint_x,PointExtrapolatedRPCFrame.y(),0.));
			LocalError RollError = rollasociated->localError(stripp);

		   std::cout<<"RPC Id: "<<rpcId<<" PointExtrapolatedRPCFrame: "<<LocalDTPoint_x<<" PointExtrapolatedRPCFrame Error: "<<LocalDTPoint_x_error<<std::endl;
  		   std::cout<<" RPC Related Local Points"<<" LocalRecHit_x: "<<LocalRecHit_x<<" LocalRecHit_x_error: "<<LocalRecHit_x_error<<" LocalRecHit_y: "<<LocalRecHit_y<<" LocalRecHit_y_error: "<<LocalRecHit_y_error<<" LocalRecHit_xy_error: "<<LocalRecHit_xy_error<<std::endl;
			

			//std::cout<<"DT RPC Id: "<<rpcId<<" Cluster Size "<<cluSize<<" Strip Width "<<stripw<<" Error: "<<cluSize*stripw/sqrt(12) <<std::endl;
			std::cout<<"DT PointExtrapolatedRPCFrameError: "<<LocalDTPoint_x_error<<" RecHitPos_Error: "<<LocalRecHit_x_error<<" Calculated error "<<error<<" RollError: "<<RollError<<std::endl;
			std::cout<<"DT Local Difference: "<<fabs(PointExtrapolatedRPCFrame.x()-recHitPos.x())<<std::endl;
			//std::cout<<"DT region  Strip Width: "<<stripw<<" RollError: "<<RollError<<" Error Calculated: "<<error<<std::endl;
			
			//Barrel
			
			
			// Station 1 Ring 0
			
		//	if( region == 0 && station ==1 && ring == 1 ){
			
					if(  layer == 1){
						//std::cout<<"Layer 1 "<<region<<" ring "<<ring<<" station "<<station<<" sector "<<sector<<" layer "<<layer<<std::endl;
						//std::cout<<"RecHit Pos "<<LocalRecHit_x<<" Error: "<<LocalRecHit_x_error<<std::endl;
						//std::cout<<"Cluster Size "<<cluSize<<" Strip Width "<<stripw<<" Error: "<<cluSize*stripw/sqrt(12) <<std::endl;
						MB1_R0_Sec1_Lay1_Local_RecHits->Fill(LocalRecHit_x);
						MB1_R0_Sec1_Lay1_Local_RecHits_Error->Fill(LocalRecHit_x_error);
						MB1_R0_Sec1_Lay1_Local_DTPoints->Fill(LocalDTPoint_x);
						MB1_R0_Sec1_Lay1_Local_DTPoints_Error->Fill(LocalDTPoint_x_error);
						MB1_R0_Sec1_Lay1_Local_Difference->Fill( Local_diff_x );
						}
						
					if( sector ==1 && layer == 2){
						//std::cout<<"Barrel Station 1 Ring 0 Sector 1 Layer 2"<<std::endl;
						//std::cout<<" RecHit Pos "<<LocalRecHit_x<<" Error: "<<LocalRecHit_x_error<<std::endl;
						//std::cout<<" Cluster Size "<<cluSize<<" Strip Width "<<stripw<<" Error: "<<error <<std::endl;
						MB1_R0_Sec1_Lay2_Local_RecHits->Fill(LocalRecHit_x);
						MB1_R0_Sec1_Lay2_Local_RecHits_Error->Fill(LocalRecHit_x_error);
						MB1_R0_Sec1_Lay2_Local_DTPoints->Fill(LocalDTPoint_x);
						MB1_R0_Sec1_Lay2_Local_DTPoints_Error->Fill(LocalDTPoint_x_error);
						MB1_R0_Sec1_Lay2_Local_Difference->Fill(Local_diff_x);
					}
					
					if( sector ==2 && layer == 1){
						MB1_R0_Sec2_Lay1_Local_RecHits->Fill(LocalRecHit_x);
						MB1_R0_Sec2_Lay1_Local_RecHits_Error->Fill(LocalRecHit_x_error);
						MB1_R0_Sec2_Lay1_Local_DTPoints->Fill(LocalDTPoint_x);
						MB1_R0_Sec2_Lay1_Local_DTPoints_Error->Fill(LocalDTPoint_x_error);
						MB1_R0_Sec2_Lay1_Local_Difference->Fill(Local_diff_x);
						}
						
					if( sector ==2 && layer == 2){
						MB1_R0_Sec2_Lay2_Local_RecHits->Fill(LocalRecHit_x);
						MB1_R0_Sec2_Lay2_Local_RecHits_Error->Fill(LocalRecHit_x_error);
						MB1_R0_Sec2_Lay2_Local_DTPoints->Fill(LocalDTPoint_x);
						MB1_R0_Sec2_Lay2_Local_DTPoints_Error->Fill(LocalDTPoint_x_error);
						MB1_R0_Sec2_Lay2_Local_Difference->Fill(Local_diff_x);
					}
				
				
				//} 			// Station 1 Ring 0
			
			}



	}		//DT Points
	}		//DT Points if
	
	
		  
// this line



//********************	CSC Points	********************	 
//********************	EndCap	********************	    	    	 

  if(CSCPoints.isValid()) if(CSCPoints->begin()!=CSCPoints->end()){//No Empty Predictions   CSC If
		RPCRecHitCollection::const_iterator rpcCSCPoint;
		
		for(rpcCSCPoint = CSCPoints->begin(); rpcCSCPoint != CSCPoints->end(); rpcCSCPoint++){	//CSC Points
				LocalPoint PointExtrapolatedRPCFrame = rpcCSCPoint->localPosition();
				LocalError PointExtrapolatedRPCFrameError = rpcCSCPoint->localPositionError();
				double LocalCSCPoint_x = PointExtrapolatedRPCFrame.x();
				double LocalCSCPoint_x_error = PointExtrapolatedRPCFrameError.xx();
				RPCDetId  rpcId = rpcCSCPoint->rpcId();
				
				std::cout<<"CSC Part: "<<std::endl;
				std::cout<<"CSCPointExtrapolatedRPCFrame: "<<LocalCSCPoint_x<<" CSC Error: "<<LocalCSCPoint_x_error<<std::endl;
 				
 				//const RPCRoll * rollasociated = rpcGeom->roll(rpcId);
      			//const TrapezoidalStripTopology* top_=dynamic_cast<const TrapezoidalStripTopology*>(&(rollasociated->topology()));
		        //LocalPoint xmin = top_->localPosition(0.);
		        
		       // std::cout<<"CSC PointExtrapolatedRPCFrame: "<<LocalCSCPoint_x<<" RecHitPos: "<<LocalRecHit_x<<std::endl;
			//std::cout<<"DT RPC Id: "<<rpcId<<" Cluster Size "<<cluSize<<" Strip Width "<<stripw<<" Error: "<<cluSize*stripw/sqrt(12) <<std::endl;
			
			
			typedef std::pair<RPCRecHitCollection::const_iterator, RPCRecHitCollection::const_iterator> rangeRecHits;
			rangeRecHits recHitCollection = rpcRecHits->get(rpcId);
			RPCRecHitCollection::const_iterator cscrecHit;	
				
				for(cscrecHit = recHitCollection.first; cscrecHit != recHitCollection.second ; cscrecHit++) {
				//countRecHits++;
				LocalPoint recHitPos=cscrecHit->localPosition();
				LocalError recHitPosError=cscrecHit->localPositionError();
				double LocalRecHit_x = recHitPos.x();
				double LocalRecHit_x_error = recHitPosError.xx();
				double LocalRecHit_y = recHitPos.y();
				double LocalRecHit_y_error = recHitPosError.yy();
				double LocalRecHit_xy_error = recHitPosError.xy();
				
				//double Local_diff_x = fabs(LocalRecHit_x - LocalCSCPoint_x );
				//int region = rpcId.region();
				//int ring = rpcId.ring();
				//int station = rpcId.station();
				//int layer = rpcId.layer();
				//int sector = rpcId.sector();
				
				const RPCRoll* rollasociated = rpcGeom->roll(rpcId);
				const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&(rollasociated->topology()));
				float stripw = top_->pitch();
				//float cluSize = cscrecHit->clusterSize();
				float error = stripw * stripw / sqrt(12);
				const float stripp = rollasociated->strip(LocalPoint(PointExtrapolatedRPCFrame.x(),PointExtrapolatedRPCFrame.y(),0.));
				LocalError RollError = rollasociated->localError(stripp);
				
				
				 std::cout<<" RPC Related Local Points"<<" LocalRecHit_x: "<<LocalRecHit_x<<" LocalRecHit_x_error: "<<LocalRecHit_x_error<<" LocalRecHit_y: "<<LocalRecHit_y<<" LocalRecHit_y_error: "<<LocalRecHit_y_error<<" LocalRecHit_xy_error: "<<LocalRecHit_xy_error<<std::endl;
			
			std::cout<<"CSC PointExtrapolatedRPCFrame: "<<LocalCSCPoint_x<<" RecHitPos: "<<LocalRecHit_x<<std::endl;
			//std::cout<<"CSC RPC Id: "<<rpcId<<" Cluster Size "<<cluSize<<" Strip Width "<<stripw<<" Error: "<<cluSize*stripw/sqrt(12) <<std::endl;
			std::cout<<"CSC PointExtrapolatedRPCFrameError: "<<LocalCSCPoint_x_error<<" RecHitPos_Error: "<<LocalRecHit_x_error<<" Calculated error "<<error<<" Roll Error "<<RollError<<std::endl;
			std::cout<<"CSC Local Difference: "<<fabs(PointExtrapolatedRPCFrame.x()-recHitPos.x())<<std::endl;
			std::cout<<"CSC region  Strip Width: "<<stripw<<" RollError: "<<RollError<<" Error Calculated: "<<error<<std::endl;
				
				}
					
			
		

		
		
		
		
		}//CSC Points
}//CSC If



				    
	
		    

}


// ------------ method called once each job just before starting event loop  ------------
void 
RPCPointAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RPCPointAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
RPCPointAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
}


// ------------ method called when ending the processing of a run  ------------
/*
void 
RPCPointAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
RPCPointAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
RPCPointAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RPCPointAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCPointAnalyzer);


/*	RPCRecHitCollection::const_iterator recHit;
    for (recHit = RecHits->begin(); recHit != RecHits->end(); recHit++) {	  
        RPCDetId rollId = (RPCDetId)(*recHit).rpcId();
        LocalPoint recHitPos=recHit->localPosition();
        const RPCRoll* rollasociated = rpcGeom->roll(rollId);
        //const BoundPlane & RPCSurface = rollasociated->surface(); 
        int region = rollId.region();
        int bx = recHit->BunchX();
		int cl = recHit->clusterSize();
		int st = recHit->firstClusterStrip();   
       // std::cout<<"RPCRecHit, RPC Id: "<<rollId<<std::endl;
        
        //	    	std::cout<<" RecHitPos "<<recHitPos<<std::endl;
        
	    if (region ==0){
 		    const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
		    float stripw = top_->pitch();
		   // double rpc_hit_pos = stripw * st ;
	    	std::cout<<"RPCRecHit Rpc Id: "<<rollId<<" Cluster Size "<<cl<<" First Cluster Strip "<<st<<" Strip Width "<<stripw<<" bx: "<<bx<<std::endl;
	    	std::cout<<" RecHitPos "<<recHitPos<<std::endl;
			   
		   }
		   
  	  }     
  */	  

