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

double beamEnergy = 2.14;
float p_mass = 0.93827;

int bin_num = 100;

float theta_min = 0;
float theta_max = 30;
float phi_min = -180;
float phi_max = 180;
float mom_min = 0;
float mom_max = 3;
float W_min = 0;
float W_max = 2;
float Q2_min = 0;
float Q2_max = 0.5;
float xB_min = 0;
float xB_max = 1.1;

double xmin = 0.9; 
double xmax = 10;

int nparts_mc = 0;
int nparts = 0;

F1D fn_sig_up = new F1D("fn_sig_up","[p0]*(x*x) + [p1]*x + [p2]", xmin, xmax);
F1D fn_sig_dn = new F1D("fn_sig_dn","[p0]*(x*x) + [p1]*x + [p2]", xmin, xmax);

// RTPC proton momenta!!
H1F h1_pmom = new H1F("h1_pmom", bin_num, 0, 0.5);
TCanvas cpmom = new TCanvas("cpmom",800,600);
cpmom.getCanvas().initTimer(1000);
cpmom.cd(0);
cpmom.draw(h1_pmom);
h1_pmom.setTitle("RTPC Proton Momentum");
h1_pmom.setTitle("momentum");

// 1D spectra
H1F h_eprime = new H1F("h_eprime", "E'", bin_num,0,2.5);
H1F h_epcal = new H1F("h_epcal", "E_{PCAL}", bin_num,0,2);
H1F h_ecin = new H1F("h_ecin", "E_{ECin}", bin_num,0,1.3);
H1F h_ecout = new H1F("h_ecout", "E_{ECout}", bin_num,0,0.5);
H1F h_ectot = new H1F("h_ectot", "E_{ECtot}", bin_num,0,4);
H1F h_nphe = new H1F("h_nphe", "nphe", 50,0,50);

// The reconstructed RGF 1D histos
H1F h1_theta = new H1F("h1_theta", bin_num, theta_min, theta_max+5);
H1F h1_phi = new H1F("h1_phi", bin_num, -phi_max, phi_max);
H1F h1_mom = new H1F("h1_mom", bin_num, 0, mom_max);
H1F h1_W = new H1F("h1_W", bin_num, W_min, W_max+0.5);
H1F h1_Q2 = new H1F("h1_Q2", bin_num, Q2_min, Q2_max);
H1F h1_xB = new H1F("h1_xB", bin_num, 0, xB_max);


H2F h2_Q2_vs_W = new H2F("h2_Q2_vs_W",100,W_min,W_max,100,Q2_min,Q2_max);
h2_Q2_vs_W.setTitle("Q2 vs W (RGF Data)");
h2_Q2_vs_W.setTitleX("W  [GeV]");
h2_Q2_vs_W.setTitleY("Q2  [GeV^2]");

H2F h2_Q2_vs_xB = new H2F("h2_Q2_vs_xB",100,xB_min,xB_max,100,Q2_min,Q2_max);
h2_Q2_vs_xB.setTitle("Q2 vs xB (RGF Data)");
h2_Q2_vs_xB.setTitleX("xB");
h2_Q2_vs_xB.setTitleY("Q2  [GeV^2]");

H2F h2_Eprime_vs_theta = new H2F("h2_Eprime_vs_theta",100,theta_min,theta_max,100,mom_min,mom_max);
h2_Eprime_vs_theta.setTitle("E' vs theta (RGF Data)");
h2_Eprime_vs_theta.setTitleX("theta");
h2_Eprime_vs_theta.setTitleY("E'  [GeV]");

H2F h2_W_vs_xB = new H2F("h2_W_vs_xB",100,xB_min,xB_max,100,W_min,W_max);
h2_W_vs_xB.setTitle("W' vs xB (RGF Data)");
h2_W_vs_xB.setTitleX("xB");
h2_W_vs_xB.setTitleY("W  [GeV]");

H2F h2_Q2_vs_W_mc = new H2F("h2_Q2_vs_W_mc",100,W_min,W_max,100,Q2_min,Q2_max);
h2_Q2_vs_W_mc.setTitle("Q2 vs W (MC Data)");
h2_Q2_vs_W_mc.setTitleX("W  [GeV]");
h2_Q2_vs_W_mc.setTitleY("Q2  [GeV^2]");

H2F h2_etot_vs_epcal = new H2F("h2_etot_vs_epcal",100,0.0,1.3,100,0.0,1.3);
h2_etot_vs_epcal.setTitle("etot vs epcal (RGF Data)");
h2_etot_vs_epcal.setTitleX("E_PCAL  [GeV]");
h2_etot_vs_epcal.setTitleY("etot [GeV]");

h1_theta.setTitleX("theta [deg]");
h1_phi.setTitleX("phi [deg]");
h1_mom.setTitleX("momentum [GeV/c]");
h1_W.setTitleX("W [GeV]");
h1_Q2.setTitleX("Q^2 [GeV^2]");
h1_xB.setTitleX("xB");

h_eprime.setTitleX("E' [GeV]");
h_epcal.setTitleX("E_PCAL [GeV]");
h_ecin.setTitleX("EC_IN [GeV]");
h_ecout.setTitleX("EC_OUT [GeV]");
h_ectot.setTitleX("EC_TOT [GeV]");
h_nphe.setTitleX("nphe");

h2_Q2_vs_W.setTitle("Q^2 vs W");
h2_Q2_vs_xB.setTitle("Q^2 vs x_B");
h2_Eprime_vs_theta.setTitle("E' vs theta");
h2_W_vs_xB.setTitle("W vs x_B");

H2F h2_sampFrac = new H2F("h2_sampFrac",200,0,9,200,0.0,0.5);
h2_sampFrac.setTitle("Sampling Fraction vs p");
h2_sampFrac.setTitleX("Momentum  [GeV]");
h2_sampFrac.setTitleY("Etot/p");

// For pre-fiducial cuts
H2F Cal_y_vs_x_precut = new H2F("Cal_y_vs_x_precut", "Cal_y_vs_x_precut", bin_num, -450,450, bin_num, -450,450);
Cal_y_vs_x_precut.setTitleX("X [cm]");
Cal_y_vs_x_precut.setTitleY("Y [cm]");

H1F Cal_lu_precut = new H1F("Cal_lu", "Cal_lu_precut", bin_num, 0, 450);
H1F Cal_lv_precut = new H1F("Cal_lv", "Cal_lv_precut", bin_num, 0, 450);
H1F Cal_lw_precut = new H1F("Cal_lw", "Cal_lw_precut", bin_num, 0, 450);

// fiducial cuts
H1F Cal_lu = new H1F("Cal_lu", "Cal_lu", bin_num, 0, 450);
H1F Cal_lv = new H1F("Cal_lv", "Cal_lv", bin_num, 0, 450);
H1F Cal_lw = new H1F("Cal_lw", "Cal_lw", bin_num, 0, 450);

H2F Cal_y_vs_x = new H2F("Cal_y_vs_x", "Cal_y_vs_x", bin_num, -450,450, bin_num, -450, 450);
Cal_y_vs_x.setTitleX("X [cm]");
Cal_y_vs_x.setTitleY("Y [cm]");

// 1D kinematic histos
TCanvas c1 = new TCanvas("c1",1000,700);
c1.getCanvas().initTimer(1000);
c1.divide(3,2);
c1.cd(0);
c1.draw(h1_theta);
c1.cd(1);
c1.draw(h1_phi);
c1.cd(2);
c1.draw(h1_mom);
c1.cd(3);
c1.draw(h1_W);
c1.cd(4);
c1.draw(h1_Q2);
c1.cd(5);
c1.draw(h1_xB);

TCanvas c2 = new TCanvas("Samp Frac & Etot",1000,600);
c2.getCanvas().initTimer(1000);
c2.divide(2,1);
c2.cd(0);
c2.draw(h2_etot_vs_epcal);
c2.cd(1);
c2.draw(h2_sampFrac);

TCanvas c3 = new TCanvas("2D Kin",1000,600);
c3.getCanvas().initTimer(1000);
c3.divide(2,1);
c3.cd(0);
c3.draw(h2_Eprime_vs_theta);
c3.cd(1);
c3.draw(h2_W_vs_xB);

TCanvas c4 = new TCanvas("2D Kin",1000,600);
c4.getCanvas().initTimer(1000);
c4.divide(2,1);
c4.cd(0);
c4.draw(h2_Q2_vs_W);
c4.cd(1);
c4.draw(h2_Q2_vs_xB);

// 1D kinematic histos
TCanvas c5 = new TCanvas("c5",1000,700);
c5.getCanvas().initTimer(1000);
c5.divide(3,2);
c5.cd(0);
c5.draw(h_eprime);
c5.cd(1);
c5.draw(h_epcal);
c5.cd(2);
c5.draw(h_ecin);
c5.cd(3);
c5.draw(h_ecout);
c5.cd(4);
c5.draw(h_ectot);
c5.cd(5);
c5.draw(h_nphe);

TCanvas can_ecal = new TCanvas("can", 1100, 600);
can_ecal.getCanvas().initTimer(1000);
can_ecal.divide(4,2);
can_ecal.cd(0);
can_ecal.draw(Cal_lu_precut);
can_ecal.cd(1);
can_ecal.draw(Cal_lv_precut);
can_ecal.cd(2);
can_ecal.draw(Cal_lw_precut);
can_ecal.cd(3);
can_ecal.draw(Cal_y_vs_x_precut);
can_ecal.cd(4);
can_ecal.draw(Cal_lu);
can_ecal.cd(5);
can_ecal.draw(Cal_lv);
can_ecal.cd(6);
can_ecal.draw(Cal_lw);
can_ecal.cd(7);
can_ecal.draw(Cal_y_vs_x);
can_ecal.save("figs/2d_ecal.png");


int cut_none = 0;
int cut_q = 0;
int cut_htcc = 0;
int cut_eprime = 0;
int cut_theta = 0;
int cut_W = 0;
int cut_Q2 = 0;
int cut_samp = 0;
int cut_fid = 0;
int cut_ec = 0;

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
    
    Bank rtpc_bank = new Bank(reader.getSchemaFactory().getSchema("RTPC::tracks"));

    Event event = new Event();
    
    while (reader.hasNext() == true) {
        parts.reset();
        calos.reset();
        chers.reset();
        
        rtpc_bank.reset();
        
        event.reset();
        
        reader.nextEvent(event);
            
        //fills bank schema with the data from this event
        event.read(parts);
        event.read(calos);
        event.read(chers);
        
        event.read(rtpc_bank);
        
        // RGA Parameters
        fn_sig_up.setParameter(0,0.0006);
        fn_sig_up.setParameter(1,-0.0058);
        fn_sig_up.setParameter(2,0.3014);
        
        fn_sig_dn.setParameter(0,-0.0018);
        fn_sig_dn.setParameter(1,0.0176);
        fn_sig_dn.setParameter(2,0.1417);
        
        PhysicsEvent physEvent = setPhysicsEvent(beamEnergy, event.read(parts));
    
        int countE = physEvent.countByPid(11); // number of electrons
        
        //System.out.println("countE: " + countE);
        
        //int max_parts = parts.getNodeLength();
        //int max_calos = calos.getNodeLength();
        //int max_chers = chers.getNodeLength();
        int max_parts = parts.getRows();
        int max_calos = calos.getRows();
        int max_chers = chers.getRows();
        
        for(int i=0; i<rtpc_bank.getRows(); i++) {
            final float momx   = rtpc_bank.getFloat("px",i);
            final float momy   = rtpc_bank.getFloat("py",i);
            final float momz   = rtpc_bank.getFloat("pz",i);
            final float pmom = Math.sqrt(momx*momx+momy*momy+momz*momz);
            
            h1_pmom.fill(pmom);
        }
        
        //System.out.println("num of particle rows: " + max_parts + ", num of calo rows: " + max_calos + ", num of cher rows: " + max_chers);
        cut_none += max_parts;
        for (int ipart=0; ipart<max_parts; ipart++) {
            //cut_none++;
    
            // require negative charge:
            if (parts.getByte("charge",ipart) < 0) { 
                cut_q++;
    
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
                
                Cal_y_vs_x_precut.fill(x_cal,y_cal);
                Cal_lu_precut.fill(lu);
                Cal_lu_precut.setLineColor(2);
                Cal_lv_precut.fill(lv);
                Cal_lv_precut.setLineColor(2);
                Cal_lw_precut.fill(lw);
                Cal_lw_precut.setLineColor(2);
                         
                if (energy>0 && nphe>5) {
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
                            
                    // ************************************************************** 
                    // ************************* begin cuts *************************
                    // **************************************************************
                               
                    if (theta > 5 && theta < 40 &&
                    e_prime > 0.1*beamEnergy && 
                    sf > fn_sig_dn.evaluate(mom) && 
                    sf <  fn_sig_up.evaluate(mom) &&
                    lu < 350 && lu > 60 && lv < 370 && lw < 390 &&
                    epcal > 0.06 && ecin > 0.025 && ecout > 0.05 && ectot > 0 ) {            
                            
                        // fill histos    
                        h_nphe.fill(nphe);
                        h_nphe.setLineColor(2);
                        
                        Cal_lu.fill(lu);
                        Cal_lu.setLineColor(2);
                        Cal_lv.fill(lv);
                        Cal_lv.setLineColor(2);
                        Cal_lw.fill(lw);
                        Cal_lw.setLineColor(2);
                        Cal_y_vs_x.fill(x_cal,y_cal);
                                        
                        h_epcal.fill(epcal);
                        h_epcal.setLineColor(2);
                        h_ecin.fill(ecin);
                        h_ecin.setLineColor(2);
                        h_ecout.fill(ecout);
                        h_ecout.setLineColor(2);
                        h_ectot.fill(ectot);
                        h_ectot.setLineColor(2);
                                    
                        h2_sampFrac.fill(mom,sf);
                        
                        h1_theta.fill(theta);
                        h1_theta.setLineColor(2);
                        h1_phi.fill(phi);
                        h1_phi.setLineColor(2);
                        h1_mom.fill(mom);
                        h1_mom.setLineColor(2);
                        h1_W.fill(W);
                        h1_W.setLineColor(2);
                        h1_Q2.fill(Q2);
                        h1_Q2.setLineColor(2);
                        h1_xB.fill(xB);
                        h1_xB.setLineColor(2);
                                    
                        h_eprime.fill(e_prime);
                        h_eprime.setLineColor(2);
                                    
                        h2_Eprime_vs_theta.fill(theta,e_prime);
                        h2_Q2_vs_W.fill(W,Q2);
                        h2_Q2_vs_xB.fill(xB,Q2);
                        h2_W_vs_xB.fill(xB,W);
                                    
                        h2_etot_vs_epcal.fill(epcal,ectot);
                                    
                        nparts++;
                    }
                                 
                } // end nphe and energy cuts
            } // end charge cut
        } // end particle loop   
    } // end while(event)
    reader.close();
} // end new line

// Normalize all the histos to 1
h1_theta.normalize(h1_theta.integral());
h1_phi.normalize(h1_phi.integral());
h1_mom.normalize(h1_mom.integral());
h1_W.normalize(h1_W.integral());
h1_Q2.normalize(h1_Q2.integral());
h1_xB.normalize(h1_xB.integral());
h_eprime.normalize(h_eprime.integral());
h_epcal.normalize(h_epcal.integral());
h_ecin.normalize(h_ecin.integral());
h_ecout.normalize(h_ecout.integral());
h_ectot.normalize(h_ectot.integral());
h_nphe.normalize(h_nphe.integral());
Cal_lu_precut.normalize(Cal_lu_precut.integral());
Cal_lv_precut.normalize(Cal_lv_precut.integral());
Cal_lw_precut.normalize(Cal_lw_precut.integral());
Cal_lu.normalize(Cal_lu.integral());
Cal_lv.normalize(Cal_lv.integral());
Cal_lw.normalize(Cal_lw.integral());


c1.save("figs/kinematics.png");
c2.save("figs/sampFrac.png");
c3.save("figs/2Dkin1.png");
c4.save("figs/2Dkin2.png");
c5.save("figs/energies.png");
can_ecal.save("figs/ecal.png");

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