// RGF e- analysis script
// 
// By: Nate Dzbenski
//
// For use, run
// rungroovy rgf_el_analysis.groovy file_list
//
// where rungroovy is aliased to where run-groovy is
// and file_list is a list of files to input

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;


import org.jlab.clas.physics.*;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.*;
import org.jlab.groot.ui.*;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.math.F1D;

GStyle.getAxisAttributesX().setTitleFontSize(32);
GStyle.getAxisAttributesY().setTitleFontSize(32);
GStyle.getAxisAttributesX().setLabelFontSize(24);
GStyle.getAxisAttributesY().setLabelFontSize(24);
GStyle.getAxisAttributesZ().setLabelFontSize(18);

double beamEnergy = 10.4;
float p_mass = 0.93827;

int bin_num = 100;

float theta_min = 0;
float theta_max = 50;
float phi_min = -180;
float phi_max = 180;
float mom_min = 0;
float mom_max = 15;
float W_min = 0;
float W_max = 15;
float Q2_min = 0;
float Q2_max = 15;
float xB_min = 0;
float xB_max = 1.1;

double xmin = 0.9; 
double xmax = 10;

int nparts_mc = 0;
int nparts = 0;

F1D fn_sig_up = new F1D("fn_sig_up","[p0]*(x*x) + [p1]*x + [p2]", xmin, xmax);
F1D fn_sig_dn = new F1D("fn_sig_dn","[p0]*(x*x) + [p1]*x + [p2]", xmin, xmax);

// Initiate histos
H1F h1_numtracks = new H1F("h1_numtracks", bin_num, -1, 6);
h1_numtracks.setTitleX("Number of tracks");

H1F h1_numhits = new H1F("h1_numhits", bin_num, -1, 100);
h1_numhits.setTitleX("Number of pads/track");

H2F h2_vze_vs_vzp = new H2F("h2_vze_vs_vzp",100,-20.0,20.0,100,-20.0,20.0);
h2_vze_vs_vzp.setTitle("Vz_e vs Vz_p");
h2_vze_vs_vzp.setTitleX("Vz_p  [cm]");
h2_vze_vs_vzp.setTitleY("Vz_e [cm]");

H2F h2_phie_vs_phip = new H2F("h2_phie_vs_phip",100,-180.0,180.0,100,-180.0,180.0);
h2_phie_vs_phip.setTitle("phi_e vs phi_p");
h2_phie_vs_phip.setTitleX("phi_p  [deg]");
h2_phie_vs_phip.setTitleY("phi_e [deg]");

H1F h1_vzdiff = new H1F("h1_vzdiff", bin_num, -20, 20);
h1_vzdiff.setTitleX("delta_vz [cm]");

H1F h1_phidiff = new H1F("h1_phidiff", bin_num, -180, 180);
h1_phidiff.setTitleX("delta_phi [deg]");

H1F h1_W = new H1F("h1_W", bin_num, 0.7, 2.0);
h1_W.setTitleX("W [GeV/c^2]");

H1F h1_Q2 = new H1F("h1_Q2", bin_num, 0.0, 0.50);
h1_Q2.setTitleX("Q^2 [GeV^2/c^2]");

H1F h1_vze = new H1F("h1_vze", bin_num, -60, 30);
h1_vze.setTitleX("vz_e [cm]");

H1F h1_theta = new H1F("h1_theta", bin_num, theta_min, theta_max);
h1_theta.setTitleX("theta [deg]");

H1F h1_phi = new H1F("h1_phi", bin_num, phi_min, phi_max);
h1_phi.setTitleX("phi [deg]");

H1F h1_emom = new H1F("h1_emom", bin_num, mom_min, mom_max);
h1_emom.setTitleX("momentum [GeV/c]");

H1F h1_xB = new H1F("h1_xB", bin_num, xB_min, xBmax);
h1_xB.setTitleX("xB");

// Initiate canvases
TCanvas c_ekin = new TCanvas("c_ekin", 1100, 600);
c_ekin.getCanvas().initTimer(1000);
c_ekin.divide(3,2);
c_ekin.cd(0);
c_ekin.draw(h1_theta);
c_ekin.cd(1);
c_ekin.draw(h1_phi);
c_ekin.cd(2);
c_ekin.draw(h1_emom);
c_ekin.cd(3);
c_ekin.draw(h1_W);
c_ekin.cd(4);
c_ekin.draw(h1_Q2);
c_ekin.cd(5);
c_ekin.draw(h1_xB);

TCanvas ctracknum = new TCanvas("ctracknum", 1100, 600);
ctracknum.getCanvas().initTimer(1000);
ctracknum.divide(2,1);
ctracknum.cd(0);
ctracknum.draw(h1_numtracks);
ctracknum.cd(1);
ctracknum.draw(h1_numhits);

TCanvas c_phi = new TCanvas("c_phi", 1100, 600);
c_phi.getCanvas().initTimer(1000);
c_phi.divide(2,1);
c_phi.cd(0);
c_phi.draw(h2_phie_vs_phip);
c_phi.cd(1);
c_phi.draw(h1_phidiff);

TCanvas c_vz = new TCanvas("c_vz", 1100, 600);
c_vz.getCanvas().initTimer(1000);
c_vz.divide(2,1);
c_vz.cd(0);
c_vz.draw(h2_vze_vs_vzp);
c_vz.cd(1);
c_vz.draw(h1_vzdiff);


// Define beam lorentz vector and target lorentz vector
LorentzVector   beam = new LorentzVector(0.0,0.0,beamEnergy,beamEnergy);
LorentzVector target = new LorentzVector(0.0,0.0,0.0,p_mass);

int lineNo = 0;

// RGF files to open
new File('.', args[0]).eachLine { line ->
    HipoReader reader=new HipoReader();
    reader.open(line);
    
    //new bank definition
    Bank parts = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
    Bank calos = new Bank(reader.getSchemaFactory().getSchema("REC::Calorimeter"));
    Bank chers = new Bank(reader.getSchemaFactory().getSchema("REC::Cherenkov"));
    
    Bank rtpc_tracks = new Bank(reader.getSchemaFactory().getSchema("RTPC::tracks"));
    Bank rtpc_hits = new Bank(reader.getSchemaFactory().getSchema("RTPC::hits"));

    Event event = new Event();
    
    while (reader.hasNext() == true) {
        parts.reset();
        calos.reset();
        chers.reset();
        
        rtpc_tracks.reset();
        rtpc_hits.reset();
        
        event.reset();
        
        reader.nextEvent(event);
            
        //fills bank schema with the data from this event
        event.read(parts);
        event.read(calos);
        event.read(chers);
        
        event.read(rtpc_tracks);
        event.read(rtpc_hits);
        
        // RGA Parameters
        fn_sig_up.setParameter(0,0.0006);
        fn_sig_up.setParameter(1,-0.0058);
        fn_sig_up.setParameter(2,0.3014);
        
        fn_sig_dn.setParameter(0,-0.0018);
        fn_sig_dn.setParameter(1,0.0176);
        fn_sig_dn.setParameter(2,0.1417);
        
        PhysicsEvent physEvent = setPhysicsEvent(beamEnergy, event.read(parts));
    
        int countE = physEvent.countByPid(11); // number of electrons
        
        int max_parts = parts.getRows();
        int max_calos = calos.getRows();
        int max_chers = chers.getRows();
        
        
        //System.out.println("num of particle rows: " + max_parts + ", num of calo rows: " + max_calos + ", num of cher rows: " + max_chers);
        
        for (int ipart=0; ipart<max_parts; ipart++) {
    
            // require negative charge:
            if (parts.getInt("pid",ipart) == 11) { 
    
                // get EC energy:
                double energy=0;
                double ecin=0;
                double epcal=0;
                double ecout=0;
                double ectot=0;
                double nphe=0;
                float x_cal = 0;
                float y_cal = 0;
                                    
                float lu = 0;
                float lv = 0;
                float lw = 0;
                    
                // get the number of photoelectrons from the cherenkov counter    
                for (int cher = 0; cher < max_chers; cher++){
                    if (calos.getShort("pindex",cher) == ipart) {
                        if (chers.getByte("detector",cher) == 15) {
                            nphe = chers.getFloat("nphe",cher);
                            break;
                        }
                     }    
                 }
                 
                 // find some calorimeter data   
                 for(int icalo = 0; icalo < max_calos; icalo++){
                     // get particle number 
                     if (calos.getShort("pindex",icalo) == ipart) {
                         if (calos.getByte("detector",icalo) == 7) {
                             if (calos.getByte("layer",icalo) == 1){
                                 epcal += calos.getFloat("energy",icalo);
                             }
                             if (calos.getByte("layer",icalo) == 4){
                                 ecin += calos.getFloat("energy",icalo);
                             }
                             if (calos.getByte("layer",icalo) == 7){
                                 ecout += calos.getFloat("energy",icalo);
                             }
                             
                             x_cal = calos.getFloat("x",icalo);
                             y_cal = calos.getFloat("y",icalo);
                                    
                             lu = calos.getFloat("lu",icalo);
                             lv = calos.getFloat("lv",icalo);
                             lw = calos.getFloat("lw",icalo);
                         }
                         
                             
                    }            
                             
                }         
                    
                energy = epcal + ecin + ecout;
                
                         
                if (energy>0) {
                
                    Particle electron = physEvent.getParticleByPid(11,0); // pid = 11, skip = 0
                            
                    final float px   = parts.getFloat("px",ipart);
                    final float py   = parts.getFloat("py",ipart);
                    final float pz   = parts.getFloat("pz",ipart);
                    final float p = Math.sqrt(px*px+py*py+pz*pz);
                            
                    // Copy lorentz vector from beam into a new LorentzVector
                    // since we are going to modify it
                    LorentzVector vecW = new LorentzVector();
                    LorentzVector vecQ2 = new LorentzVector();
                            
                    // creates a copy of lorentz vector from electron
                    LorentzVector  vecE = electron.vector();
                    
                    vecW.copy(target);
                    vecW.add(vecE);
                    vecQ2.copy(beam);
                    vecQ2.sub(vecE);
                            
                    double e_prime = vecE.e();
                    double theta = vecE.theta();
                    double phi = vecE.phi();
                            
                    theta *= 180/Math.PI;
                    phi *= 180/Math.PI;
                    
                    double nu = beamEnergy - e_prime;
                            
                    double mom = vecE.p();
                    double W = Math.sqrt(p_mass*p_mass + vecQ2.mass2() + 2*p_mass*nu);
                    double Q2 = -vecQ2.mass2();
                    double xB = Q2/(2.0*p_mass*nu);
                            
                    ectot = ecin+ecout;
                             
                    double sf = energy/mom;
                    
                    double e_vz = parts.getFloat("vz",ipart);
                    
                    // fill electron kinematic histos 
                    h1_W.fill(W);
                    h1_Q2.fill(Q2);
                    h1_vze.fill(e_vz);
                    h1_xB.fill(xB);
                    h1_theta.fill(theta);
                    h1_phi.fill(phi);
                    h1_emom.fill(p);
                    
                    if (!rtpc_hits) {
                        System.out.println("No RTPC hits!")
                        continue;
                    } 
            
                    // ************************************************************** 
                    // ************************* begin cuts *************************
                    // **************************************************************
                    //if(event.hasBank("RTPC::hits")){System.out.println("We have a bank!!");}           
                    //if (Q2 > 0.05 && Q2 < 0.1 && W > 0.85 && W < 1.05 && e_vz > -20.0 && e_vz < 15.0) {            
                        
                        // Let's look into the rtpc bank
                        int num_rtpc_tracks = rtpc_tracks.getRows();
                        int num_rtpc_hits = rtpc_hits.getRows();
                        
                        int tid = 0;
                        int _tid = -991;
                        int pads_per_track = 0;
                        
                            
                            h1_numtracks.fill(num_rtpc_tracks);
                            
                            for(int k = 0; k < num_rtpc_hits; k++){
                                if(k == 0){ 
                                    _tid = rtpc_hits.getInt("trkID",k);
                                    pads_per_track++;
                                }
                                else {
                                    tid = rtpc_hits.getInt("trkID",k);
                                    if(tid == _tid) {pads_per_track++;}
                                    else {
                                        h1_numhits.fill(pads_per_track);
                                        _tid = tid;
                                        pads_per_track = 0;
                                    }
                                }
                            }
                            
                        if(num_rtpc_tracks > 0 && num_rtpc_hits > 0){
                            for(int itr = 0; itr < num_rtpc_tracks; itr++){
                                float momx   = rtpc_tracks.getFloat("px",itr);
                                float momy   = rtpc_tracks.getFloat("py",itr);
                                float momz   = rtpc_tracks.getFloat("pz",itr);
                                float pmom = Math.sqrt(momx*momx+momy*momy+momz*momz);
                                
                                double p_phi = Math.atan2(momy,momx);
                                double e_phi = vecE.phi();
                                e_phi *=  180/Math.PI;
                                p_phi *=  180/Math.PI;
                                
                                double p_vz = rtpc_tracks.getFloat("vz",itr);
                                
                                h2_vze_vs_vzp.fill(p_vz, e_vz);
                                h2_phie_vs_phip.fill(p_phi, e_phi);
    
                                h1_vzdiff.fill(e_vz-p_vz);
                                h1_phidiff.fill(e_phi-p_phi);
                                
                            }
                           
                            
                        }
                         
                    //}
                                 
                } // end nphe and energy cuts
            } // end charge cut
        } // end particle loop   
    } // end while(event)
    reader.close();
} // end new line


ctracknum.save("figs/proton/track_info.png");
c_ekin.save("figs/electron/ekinematics.png");
c_phi.save("figs/proton/phi.png");
c_vz.save("figs/proton/vz.png");

// defining method because getPhysicsEvent only works for one type of bank
public static PhysicsEvent setPhysicsEvent(double beam, Bank parts) {

    PhysicsEvent physEvent = new PhysicsEvent();
    physEvent.setBeamParticle(new Particle(11, 0.0D, 0.0D, beam));
 
    if (!parts) {
        System.out.println("No event!")
        return physEvent;
    } 
    
    else {
        for (int i = 0; i < parts.getRows(); ++i) {
            final float px   = parts.getFloat("px",i);
            final float py   = parts.getFloat("py",i);
            final float pz   = parts.getFloat("pz",i);
            final float vx   = parts.getFloat("vx",i);
            final float vy   = parts.getFloat("vy",i);
            final float vz   = parts.getFloat("vz",i);
                        
            int pid = parts.getInt("pid",i);
            
            int status = parts.getInt("status",i);
            int detector = 1;
            if (status >= 2000 && status < 3000) {
                detector = 2;
            }

            if (status >= 4000) {
                detector = 3;
            }

            Particle p = new Particle();
            if (pid != 0) {
                p.initParticle(pid, (double) px, (double) py, (double) pz, (double) vx, (double) vy, (double) vz);
            } 
            else {
                p.initParticleWithPidMassSquare(pid, 0, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D);
            }

            //p.setStatus(detector);
            physEvent.addParticle(p);
        }

        return physEvent;
    }
}