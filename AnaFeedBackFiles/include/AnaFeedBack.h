#include "EdbSegP.h"
#include "EdbMomentumEstimator.h"
#include "EdbPattern.h"
#include "EdbSegCouple.h"
#include "EdbPVRec.h"
#include "EdbVertex.h"
#include "EdbDataSet.h"


void ParseFeedBack(std::string filename){

    char FlagPart[][10]={"", "mu", "charm", "e", "e-pair"};
	char FlagCS[][10]   ={"", "cs"};

    FILE *fp = fopen(filename.c_str(),"rt");
	if(fp==NULL){
		printf(Form("File open error : %s . stop.\n", filename));
		return;
	}
	
	EdbPVRec *pvr = new EdbPVRec();
	EdbScanCond cond;
	EdbTrackP *t=0;

    char buf[256];

    while(fgets(buf,sizeof(buf),fp)){
		TString line=buf;
		printf("%s", line.Data());
		
		// Remove comments
		line.ReplaceAll("\t"," ");
		int iposcmt = line.Index("//");
		if(iposcmt!=kNPOS) line.Remove(iposcmt, line.Length()-iposcmt);
		
		// Check number of columns.
		TObjArray *tokens = line.Tokenize(" ");
		int ntokens = tokens->GetEntries();
		delete tokens;
		
		
		if( ntokens==0 ) continue;
		
		else if( ntokens == 10 ){
			// Vertices
			float x,y,z; int id,isprimary,istau, ischarm;
			sscanf(line.Data(),"%d %f %f %f %d %d %d", &id, &x, &y, &z, &isprimary, &ischarm, &istau);
			EdbVertex *v = new EdbVertex();
			v->SetXYZ(x,y,z); 
			v->SetID(id);
			v->SetFlag(isprimary);
			pvr->AddVertex(v);
			printf("Vertex %d %f %f %f\n",v->ID(), v->X(), v->Y(), v->Z());
		}

		else if( ntokens == 27 ){
			// Tracks
			float x,y,z,ax,ay, ip_upvtx, ip_downvtx,  p,pmin,pmax;
			int id_track, id_upvtx, id_downvtx,  manual, particle_id, scanback, darkness;
			int OutOfBrick, LastPlate, nseg, iplRmax1, iplRmax2, result;
			float RmaxT, RmaxL, rmst, rmsl;
			
			int n,nn;
			sscanf(line.Data(),         "%d %d %d %f %f %f %f %f %f%n", &id_track, &id_upvtx, &id_downvtx, &x, &y, &z, &ax, &ay, &ip_upvtx, &nn);
			sscanf(line.Data()+(n=nn),  "%f %f %f %f %d %d %d %d %d%n", &ip_downvtx, &p,&pmin,&pmax, &manual, &particle_id, &scanback, &darkness, &OutOfBrick, &nn);
			sscanf(line.Data()+(n+=nn), "%d %d %f %f %f %f %d %d %d",   &LastPlate, &nseg, &RmaxT, &RmaxL, &rmst, &rmsl, &iplRmax1, &iplRmax2, &result);
			
			// Create Track. "t" is defined out of main loop.
			t = new EdbTrackP;
			t->Set(id_track, x, y, ax, ay, 0, 0);
			t->SetZ(z);
			t->SetTrack(id_track);
			pvr->AddTrack(t);

			// fill COV for vertexing
			t->SetErrors();
			cond.FillErrorsCov(t->TX(), t->TY(), t->COV());
			
			// associate to vertex.
			for(int i=0; i<pvr->Nvtx(); i++){
				EdbVertex *v = pvr->GetVertex(i);
				if(id_upvtx==v->ID()||id_downvtx==v->ID()){
					EdbVTA *vta = new EdbVTA(t, v);
					vta->SetZpos( t->Z()>v->Z() ? 1 : 0);
					vta->SetFlag(2);
					vta->SetImp( id_upvtx==v->ID()? ip_upvtx : ip_downvtx);
					v->AddVTA(vta);
				}
			}
		}
		
		else if( ntokens ==  9 ){
			// Segments
			int ipl, type, irec, ph;
			float x,y,z,ax,ay;
			sscanf(line.Data(),"%d %f %f %f %f %f %d %d %d", &ipl, &x, &y, &z, &ax, &ay, &type, &irec, &ph);
			
			EdbSegP *s = new EdbSegP(t->ID(),x,y,ax,ay,0,0);
			s->SetZ(z);
			s->SetPID(ipl);
			s->SetPlate(ipl);
			s->SetW(ph);
			s->SetTrack(t->ID());
			s->SetSide(type);

			// fill COV for vertexing
			s->SetErrors();
			cond.FillErrorsCov(s->TX(), s->TY(), s->COV());
			
			// Add segment in PVRec and Track, keeping consistency of pointer in them.
			EdbSegP *s_in_pattern = pvr->AddSegment(*s);
			t->AddSegment(s_in_pattern);
		}
		else{
			std::cout << "Warning: Reading feedback file: items = " << ntokens << ". This is NOT specified!" << std::endl;
		}
		
	} // The end of while loop of reading feedback files line by line

    for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetPID(i);
	for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetSegmentsPID();
	for(int i=0;i<pvr->Ntracks(); i++) pvr->GetTrack(i)->SetCounters();

    printf("\nEdbPVRec summary. %d vertices, %d tracks.\n", pvr->Nvtx(), pvr->Ntracks());
	pvr->Print();

	EdbDataProc *dproc = new EdbDataProc();
	dproc->MakeTracksTree(pvr, "linked_tracks.root");
	
	fclose(fp);

}


