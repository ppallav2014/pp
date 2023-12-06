// Detector Construction class.
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
// Action initialization class.
#include "G4VUserActionInitialization.hh"
// Run Action class.
#include "G4AccumulableManager.hh"
#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4UnitsTable.hh"
#include "G4Run.hh"
// Event action class.
#include "G4UserEventAction.hh"
#include "G4Event.hh"
// Stepping action class.
#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
// Primary generator action class.
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
// Physics list.
#include "FTFP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "G4VUserPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
//#include "G4MuIonisation.hh"
//#include "G4DecayPhysics.hh"
#include "Shielding.hh"
#include "QGSP_BIC_HP.hh"
// Main.
//#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "globals.hh"
//#include "EventAction.hh"
#include <fstream>
#include <vector>
#include "G4Threading.hh"

//*****************************************************************Parameters class****************************************************//

class Parameters {
public:
    G4int numberOfEvent     = 10000;
    G4bool GUI              = false;
    G4bool checkOverlaps    = true;
    G4int NumberOfThreads   = 8;
    // All units in cm
    // Detector Size
    G4double  Det_sizeXY    = 15.0;
    G4double  Det_sizeZ     = 1;
    //Detector Position
    G4double RPCPosX[6]     = {0, 0, 2.5, 0.0, 0.0};
    G4double RPCPosY[6]     = {0, 0.0,0.0, 0.0, 0.0, 0.0};
    G4double RPCPosZ[6]     = {-6.5, -8.5, 0, 0, 0, 0};
    // Concrete Size and pos
    G4double ConcreteSize[3] = {100.0,  100.0, 100.0};
    G4double ConcretePos[3]  = { 0.0,  0.0,  0.0};
G4double BlockY = 18;
G4double BlockX = 18;
G4double BlockZ = 3.0;

std::string filename    = "../output/0.txt";
G4double BlockZpos = 8.0*0.5*cm;
//G4double BlockZpos = 16.0*0.5*cm;
      
G4double BlockXpos = 0;
G4double BlockYpos = 0;




//******************wood parameters**********
G4double wBlockY = 20;
G4double wBlockX = 20;
G4double wBlockZ = 2;


G4double wBlockXpos = 2;
G4double wBlockYpos = 0;
G4double wBlockZpos = -1*cm;
    // G4String particleName   = "mu-";
		};

//****************************************************DetectorConstruction class*****************************************************//
class DetectorConstruction : public G4VUserDetectorConstruction {
private:
    Parameters parameters;
    G4LogicalVolume*  fScoringVolume;
public:
    DetectorConstruction(): G4VUserDetectorConstruction(), fScoringVolume(0) { }
    ~DetectorConstruction(){ }
    G4VPhysicalVolume* Construct(){
        G4NistManager* nist     = G4NistManager::Instance();
        G4Material* Air         = nist->FindOrBuildMaterial("G4_AIR");
        G4Material* Argon       = nist->FindOrBuildMaterial("G4_Ar");
        G4Material* Lead        = nist->FindOrBuildMaterial("G4_Pb");
        G4Material* Uranium     = nist->FindOrBuildMaterial("G4_U");
        G4Material* Iron        = nist->FindOrBuildMaterial("G4_Fe");
        G4Material* Aluminum    = nist->FindOrBuildMaterial("G4_Al");
        G4Material* Zinc        = nist->FindOrBuildMaterial("G4_Zn");
        G4Material* Copper      = nist->FindOrBuildMaterial("G4_Cu");
        G4Material* WoodMaterial = nist->FindOrBuildMaterial("G4_W");
       
       
        G4double  Det_sizeXY    = parameters.Det_sizeXY *cm;
        G4double  Det_sizeZ     = parameters.Det_sizeZ *cm;
        // World
        G4double world_sizeXY   = (Det_sizeXY + 80.0) *cm; // must be greater than Det_sizeXY
        G4double world_sizeZ    = ((parameters.RPCPosZ[0] - parameters.RPCPosZ[5])+80.0)*cm;
        G4bool checkOverlaps    = parameters.checkOverlaps;
        // Geometry
        
	//Rotation matrix//
	G4double rotationAngle = 0.0 * degree; //for th second plate (origine)
	G4RotationMatrix* rotation = new G4RotationMatrix();
	rotation->rotateY(rotationAngle);

	// G4Box("Name", 0.5*Size_X, 0.5*Size_Y, 0.5*Size_Z);
        G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
        // G4LogicalVolume(solid_poiner, its material, "Name");
        G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World");
        // G4PVPlacement(rotation, position, logical vol, name, mother vol, boolean op, copy no, Overlaps)
        G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);
        



        // RPCilators
        G4ThreeVector Pos_1, Pos_2, Pos_3;
    	Pos_1 = G4ThreeVector(parameters.RPCPosX[0]*cm, parameters.RPCPosY[0]*cm,  parameters.RPCPosZ[0]*cm) ;//(a: top detector)
   	Pos_2 = G4ThreeVector(parameters.RPCPosX[1]*cm, parameters.RPCPosY[1]*cm,  parameters.RPCPosZ[1]*cm); //(b: middle detector).
       // Pos_3 = G4ThreeVector(parameters.RPCPosX[2]*cm, parameters.RPCPosY[2]*cm,  parameters.RPCPosZ[2]*cm);//(c: bottom detector)
        
        G4Box* solidRPC = new G4Box("RPC", Det_sizeXY*0.5, Det_sizeXY*0.5, Det_sizeZ*0.5);
        G4LogicalVolume* logicRPC = new G4LogicalVolume(solidRPC, Argon, "RPC");

     	 new G4PVPlacement(rotation, Pos_1, logicRPC, "RPC", logicWorld,  false, 0, checkOverlaps);     // RPC 1
   	 new G4PVPlacement(rotation, Pos_2, logicRPC, "RPC", logicWorld,  false, 1, checkOverlaps);     // RPC 2
        // new G4PVPlacement(rotation, Pos_3, logicRPC, "RPC", logicWorld,  false, 2, checkOverlaps);     // RPC 3

/*
   	// Lead block
	G4VSolid* solidBlock1 = new G4Box("Block1", parameters.BlockX*cm*0.5, parameters.BlockY*cm*0.5, parameters.BlockZ*cm*0.5);
	G4LogicalVolume* logicBlock1 = new G4LogicalVolume(solidBlock1, Lead, "Lead"); // Use the Lead material	
	new G4PVPlacement(0, G4ThreeVector(parameters.BlockXpos, parameters.BlockYpos, parameters.BlockZpos), logicBlock1, "Lead", logicWorld, false, 0, checkOverlaps);

// Wood block
	G4VSolid* solidBlock2 = new G4Box("Block2", parameters.wBlockX*cm*0.5, parameters.wBlockY*cm*0.5, parameters.wBlockZ*cm*0.5);
	G4LogicalVolume* logicBlock2 = new G4LogicalVolume(solidBlock2, WoodMaterial, "WoodMaterial"); // Use the Zinc material
	new G4PVPlacement(0, G4ThreeVector(parameters.wBlockXpos, parameters.wBlockYpos, parameters.wBlockZpos), logicBlock2, "WoodMaterial", logicWorld, false, 0, checkOverlaps);

*/


        // Set ScoringVolume
        fScoringVolume = logicRPC;
        //always return the physical World
        return physWorld;
    }
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
};

//**********************************************************PrimaryGeneratorAction class***********************************************************//

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
    Parameters parameters;
    G4ParticleGun*  fParticleGun;
    std::ifstream infile;
public:
    PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
        fParticleGun  = new G4ParticleGun(1);   // G4ParticleGun(G4int number_particle_each_event)
        infile.open("../CRY_Data_File.txt");  // CRY data file
    }
    ~PrimaryGeneratorAction() { delete fParticleGun; }
    void GeneratePrimaries(G4Event* anEvent) {
    G4double x0 = (parameters.BlockX )*1 *  (G4UniformRand() - 0.5);
    G4double y0 = (parameters.BlockY )*1 * (G4UniformRand() - 0.5);
    G4double z0 =  (parameters.BlockZ ) + 0.5 ; // z position -65.0 cm

    // Calculate new position after rotation
    G4double rotationAngleDegrees = 0.0; // Replace with your desired rotation angle in degrees
    G4double rotationAngleRadians = rotationAngleDegrees * M_PI / 180.0;

    // Calculate new x and z coordinates after rotation
    G4double rotatedX = x0 * cos(rotationAngleRadians) + z0 * sin(rotationAngleRadians);
    G4double rotatedZ = -x0 * sin(rotationAngleRadians) + z0 * cos(rotationAngleRadians);

    G4int PDGid;
    G4double En, xx, yy, px, py, pz;
    infile >> PDGid >> En >> xx >> yy >> px >> py >> pz;

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(PDGid);

    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(En * MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
    fParticleGun->SetParticlePosition(G4ThreeVector(rotatedX * cm, y0 * cm, rotatedZ * cm));

    fParticleGun->GeneratePrimaryVertex(anEvent);
}

};

//********************************************************************RunAction class***********************************************//

class RunAction : public G4UserRunAction {
public:
    RunAction(): G4UserRunAction(){ }
    ~RunAction(){ }
    void BeginOfRunAction(const G4Run*){ }
    void EndOfRunAction(const G4Run* run){
        G4int nofEvents = run->GetNumberOfEvent();
        if (nofEvents == 0) return;
        if (IsMaster()) G4cout << G4endl << "------> " << ".------------------.End of Global Run.----------------------." << G4endl;
        else G4cout << G4endl << "--------------------End of Local Run------------------------" << G4endl;
    }
};

//***************************************************************EventAction class*****************************//

class EventAction : public G4UserEventAction {
public:
    double RPCHits[6][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    G4double SetKE = 0.0;
    RunAction* fRunAction;
    Parameters parameters;
    std::fstream outfile;
public:
    EventAction(RunAction* runAction) : G4UserEventAction(), fRunAction(runAction) { }
    ~EventAction() { }
    void BeginOfEventAction(const G4Event*) { }
    void EndOfEventAction(const G4Event* event) {
        G4int eventID = event->GetEventID();
        int per = 1000;
        if (eventID%per == 0) {
		G4cout << "Total Events : " << eventID<< G4endl;
		}

	 outfile.open(parameters.filename, std::ios::out | std::ios::app);
	 if (RPCHits[0][0] != 0.0 || RPCHits[1][0] != 0.0 || RPCHits[2][0] != 0.0) {
        outfile << (RPCHits[0][0] != 0.0 ? 1 : 0) << " "
                << (RPCHits[1][0] != 0.0 ? 1 : 0) << " "
                << (RPCHits[2][0] != 0.0 ? 1 : 0) << G4endl;
    } else {
        outfile << "0 0 0" << G4endl;
    }
    outfile.close();

    eventID = 0;

    for (G4int i = 0; i < 6; i++) {
        for (G4int j = 0; j < 3; j++) {
            outfile << RPCHits[i][j] << " ";
        }
    }

    outfile << G4endl;
    outfile.close();

    for (G4int i = 0; i < 6; i++) {
        for (G4int j = 0; j < 3; j++) {
            RPCHits[i][j] = 0.0;
        }
    }
}

      
};

//******************************************************SteppingAction class************************************************//

class SteppingAction : public G4UserSteppingAction {
private:
    EventAction* fEventAction;
    G4LogicalVolume* fScoringVolume;
public:
    SteppingAction(EventAction* eventAction) : G4UserSteppingAction(), fEventAction(eventAction), fScoringVolume(0) { }
    ~SteppingAction() { }
    void UserSteppingAction(const G4Step* step) {
        if (!fScoringVolume) {
            const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            fScoringVolume = detectorConstruction->GetScoringVolume();
        }

        G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        if (volume == fScoringVolume) {
            G4int copy_no = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
            const G4Track* track = step->GetTrack();
            G4int pid = track->GetParticleDefinition()->GetPDGEncoding();

            if (pid == 13 || pid == -13 ) {
            
            //|| pid == 11 || pid == -11 || pid == 22 || pid == -22
            
            
            //            if (pid == 13 || pid == -13 || pid == 11 || pid == -11 || pid == 2212 || pid -2212 || pid == 22 || pid == -22) {
                G4double energy = track->GetDynamicParticle()->GetKineticEnergy();
                if (energy > 1.0   and energy < 5000) 
                fEventAction->RPCHits[copy_no][0] += energy;
                fEventAction->RPCHits[copy_no][1] += 1.0;
                fEventAction->RPCHits[copy_no][2] = 1.0;
            }
        }
    }
};

//*********************************************************************ActionInitialization class****************************************//

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(): G4VUserActionInitialization() { }
    ~ActionInitialization() { }
    void BuildForMaster() const {
        RunAction* runAction = new RunAction();
        SetUserAction(runAction);
    }
    void Build() const {
        RunAction* runAction = new RunAction();
        SetUserAction(runAction);
        EventAction* eventAction = new EventAction(runAction);
        SetUserAction(eventAction);
        SteppingAction* steppingAction = new SteppingAction(eventAction);
        SetUserAction(steppingAction);
        PrimaryGeneratorAction* primaryGeneratorAction = new PrimaryGeneratorAction();
        SetUserAction(primaryGeneratorAction);
  };  
};

//*******************************************************************Physics list*****************************************************************//

/*
class MyPhysicsConstructor : public QGSP_BERT_HP {
public:
MyPhysicsConstructor() {
RegisterPhysics(new G4EmStandardPhysics_option3()); // Electromagnetic physics
   //RegisterPhysics(new G4DecayPhysics()); // Electromagnetic physics
           // RegisterPhysics(new G4RadioactiveDecayPhysics()); // Electromagnetic physics
  }

 ~MyPhysicsConstructor() {}

};
*/
//********************************************************************************main*******************************************/


int main(int argc, char** argv) {

   Parameters parameters;
   G4bool GUI = parameters.GUI;
   G4Random::setTheEngine(new CLHEP::RanecuEngine);
   //auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
   G4RunManager* runManager = new G4RunManager;
    Shielding* Physics = new Shielding();
   runManager->SetNumberOfThreads(parameters.NumberOfThreads);
   DetectorConstruction* Detector = new DetectorConstruction();
 // MyPhysicsConstructor* Physics = new MyPhysicsConstructor();
   ActionInitialization* Action = new ActionInitialization();
   runManager->SetUserInitialization(Detector);
   runManager->SetUserInitialization(Physics);
   runManager->SetUserInitialization(Action);

   if ( ! GUI ) {
       // Start a Run
       runManager->Initialize();
       int numberOfEvent = parameters.numberOfEvent;
       runManager->BeamOn(numberOfEvent);
       delete runManager;
   } else {
       // Alternative Run
       G4UIExecutive* ui = 0;
       if ( argc == 1 ) ui = new G4UIExecutive(argc, argv);
       G4VisManager* visManager = new G4VisExecutive;
       visManager->Initialize();
       G4UImanager* UImanager = G4UImanager::GetUIpointer();
       if ( ! ui ) {
           G4String command = "/control/execute ";
           G4String fileName = argv[1];
           UImanager->ApplyCommand(command+fileName);
       }
       else {
           UImanager->ApplyCommand("/control/execute init_vis.mac");
           ui->SessionStart();
           delete ui;
       }
       delete visManager;
       delete runManager;
   }
}
