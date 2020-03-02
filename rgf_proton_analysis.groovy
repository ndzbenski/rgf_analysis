// RGF proton analysis script
// 
// By: Nate Dzbenski
//
// For use, run
// rungroovy rgf_proton_10p4.groovy file_list
//
// where rungroovy is aliased to where run-groovy is
// and file_list is a list of files to input

import javax.swing.JFrame;
import javax.swing.JTabbedPane;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

import org.jlab.clas.physics.*;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.*;
import org.jlab.groot.ui.*;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.math.F1D;
import org.jlab.groot.graphics.EmbeddedCanvas;

GStyle.getAxisAttributesX().setTitleFontSize(32);
GStyle.getAxisAttributesY().setTitleFontSize(32);
GStyle.getAxisAttributesX().setLabelFontSize(24);
GStyle.getAxisAttributesY().setLabelFontSize(24);
GStyle.getAxisAttributesZ().setLabelFontSize(18);

JFrame frame = new JFrame("RGF Analysis");
frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
frame.setSize(1200, 600);
JTabbedPane tabbedPane = new JTabbedPane(); 

JFrame pframe = new JFrame("Elastic Proton Analysis");
pframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
pframe.setSize(1200, 600);
JTabbedPane tabbedPaneP = new JTabbedPane();       

double beamEnergy = 10.4;
float p_mass = 0.93827;

int bin_num = 100;

float theta_min = 4;
float theta_max = 40;
float phi_min = -180;
float phi_max = 180;
float mom_min = 0;
float mom_max = 10;
float W_min = 0.7;
float W_max = 5.0;
float Q2_min = 0.0;
float Q2_max = 10.0;
float xB_min = 0;
float xB_max = 1.2;

double xmin = 0.9; 
double xmax = 10;

int nparts_mc = 0;
int nparts = 0;

F1D fn_sig_up = new F1D("fn_sig_up","[p0]*(x*x) + [p1]*x + [p2]", xmin, xmax);
F1D fn_sig_dn = new F1D("fn_sig_dn","[p0]*(x*x) + [p1]*x + [p2]", xmin, xmax);

// Initiate histos
H1F h1_numtracks = new H1F("h1_numtracks", bin_num, -1, 15);
H1F h1_numhits = new H1F("h1_numhits", bin_num, -1, 100);
H1F h1_momdiff = new H1F("h1_momdiff", bin_num, -0.50, 0.5);
H1F h1_thdiff = new H1F("h1_thdiff", bin_num, -90, 90);
H1F h1_vzdiff = new H1F("h1_vzdiff", bin_num, -60, 60);
H1F h1_phidiff = new H1F("h1_phidiff", bin_num, -360, 360);
H1F h1_W = new H1F("h1_W", 2*bin_num, W_min, W_max);
H1F h1_Q2 = new H1F("h1_Q2", 4*bin_num, Q2_min, Q2_max);
H1F h1_vze = new H1F("h1_vze", bin_num, -60, 30);
H1F h1_theta = new H1F("h1_theta", bin_num, theta_min, theta_max);
H1F h1_phi = new H1F("h1_phi", bin_num, phi_min, phi_max);
H1F h1_emom = new H1F("h1_emom", 2*bin_num, mom_min, mom_max);
H1F h1_xB = new H1F("h1_xB", bin_num, xB_min, xB_max);
H1F h1_tshift = new H1F("h1_tshift", bin_num, -5000, 5000);
H1F h1_pmom = new H1F("h1_pmom", bin_num, 0.0, 2.0);
H1F h1_ptheta = new H1F("h1_ptheta", bin_num, -5, 180);


H1F h1_numtracksu = new H1F("h1_numtracks_uncut", bin_num, -1, 15);
h1_numtracksu.setTitleX("Number of tracks");
h1_numtracksu.setLineColor(2);
H1F h1_numhitsu = new H1F("h1_numhits_uncut", bin_num, -1, 100);
h1_numhitsu.setTitleX("Number of pads/track");
h1_numhitsu.setLineColor(2);
H1F h1_momdiffu = new H1F("h1_momdiff_uncut", bin_num, -0.50, 0.5);
h1_momdiffu.setTitleX("delta_mom [GeV/c]");
h1_momdiffu.setLineColor(2);
H1F h1_thdiffu = new H1F("h1_thdiff_uncut", bin_num, -90, 90);
h1_thdiffu.setTitleX("delta_theta [deg]");
h1_thdiffu.setLineColor(2);
H1F h1_vzdiffu = new H1F("h1_vzdiff_uncut", bin_num, -60, 60);
h1_vzdiffu.setTitleX("delta_vz [cm]");
h1_vzdiffu.setLineColor(2);
H1F h1_phidiffu = new H1F("h1_phidiff_uncut", bin_num, -360, 360);
h1_phidiffu.setTitleX("delta_phi [deg]");
h1_phidiffu.setLineColor(2);
H1F h1_Wu = new H1F("h1_W_uncut", 2*bin_num, W_min, W_max);
h1_Wu.setTitleX("W [GeV/c^2]");
h1_Wu.setLineColor(2);
H1F h1_Q2u = new H1F("h1_Q2_uncut", 4*bin_num, Q2_min, Q2_max);
h1_Q2u.setTitleX("Q^2 [GeV^2/c^2]");
h1_Q2u.setLineColor(2);
H1F h1_vzeu = new H1F("h1_vze_uncut", bin_num, -60, 30);
h1_vzeu.setTitleX("vz_e [cm]");
h1_vzeu.setLineColor(2);
H1F h1_thetau = new H1F("h1_theta_uncut", bin_num, theta_min, theta_max);
h1_thetau.setTitleX("theta [deg]");
h1_thetau.setLineColor(2);
H1F h1_phiu = new H1F("h1_phi_uncut", bin_num, phi_min, phi_max);
h1_phiu.setTitleX("phi [deg]");
h1_phiu.setLineColor(2);
H1F h1_emomu = new H1F("h1_emom_uncut", 2*bin_num, mom_min, mom_max);
h1_emomu.setTitleX("momentum [GeV/c]");
h1_emomu.setLineColor(2);
H1F h1_xBu = new H1F("h1_xB_uncut", bin_num, xB_min, xB_max);
h1_xBu.setTitleX("xB");
h1_xBu.setLineColor(2);
H1F h1_tshiftu = new H1F("h1_tshift_uncut", bin_num, -5000, 5000);
h1_tshiftu.setTitleX("tshift [ns]");
h1_tshiftu.setLineColor(2);
H1F h1_pmomu = new H1F("h1_pmom_uncut", bin_num, 0.0, 2.0);
h1_pmomu.setTitleX("mom_p [GeV/c]");
h1_pmomu.setLineColor(2);
H1F h1_pthetau = new H1F("h1_ptheta_uncut", bin_num, -5, 180);
h1_pthetau.setTitleX("theta_p [deg]");
h1_pthetau.setLineColor(2);

H2F h2_vze_vs_vzp = new H2F("h2_vze_vs_vzp",100,-25.0,25.0,100,-30.0,22.0);
h2_vze_vs_vzp.setTitle("Vz_e vs Vz_p");
h2_vze_vs_vzp.setTitleX("Vz_p  [cm]");
h2_vze_vs_vzp.setTitleY("Vz_e [cm]");

H2F h2_phie_vs_phip = new H2F("h2_phie_vs_phip",100,-180.0,180.0,100,-180.0,180.0);
h2_phie_vs_phip.setTitle("phi_e vs phi_p");
h2_phie_vs_phip.setTitleX("phi_p  [deg]");
h2_phie_vs_phip.setTitleY("phi_e [deg]");

// For 2.14 GeV analysis
H2F h2_ptheta = new H2F("h2_ptheta",bin_num, 75,81 ,bin_num,0.0,181);
h2_ptheta.setTitle("theta_pred vs theta_meas");
h2_ptheta.setTitleX("theta_predicted  [deg]");
h2_ptheta.setTitleY("theta_measured [deg]");

H2F h2_mom = new H2F("h2_mom",bin_num,0.22,0.34,bin_num,0.0,0.50);
h2_mom.setTitle("mom_pred vs mom_meas");
h2_mom.setTitleX("mom_predicted [GeV/c]");
h2_mom.setTitleY("mom_measured [GeV/c]");

// Initiate canvases
EmbeddedCanvas c_ekin = new EmbeddedCanvas();
c_ekin.initTimer(1000);
c_ekin.divide(3,2);
c_ekin.cd(0);
c_ekin.draw(h1_thetau);
c_ekin.draw(h1_theta,"same");
c_ekin.cd(1);
c_ekin.draw(h1_phiu);
c_ekin.draw(h1_phi,"same");
c_ekin.cd(2);
c_ekin.draw(h1_emomu);
c_ekin.draw(h1_emom,"same");
c_ekin.cd(3);
c_ekin.draw(h1_Wu);
c_ekin.draw(h1_W,"same");
c_ekin.cd(4);
c_ekin.draw(h1_Q2u);
c_ekin.draw(h1_Q2,"same");
c_ekin.cd(5);
c_ekin.draw(h1_xBu);
c_ekin.draw(h1_xB,"same");

EmbeddedCanvas c_p1d = new EmbeddedCanvas();
c_p1d.initTimer(1000);
c_p1d.divide(3,1);
c_p1d.cd(0);
c_p1d.draw(h1_tshift,"same");
c_p1d.draw(h1_tshiftu);
c_p1d.cd(1);
c_p1d.draw(h1_pmomu);
c_p1d.draw(h1_pmom,"same");
c_p1d.cd(2);
c_p1d.draw(h1_pthetau);
c_p1d.draw(h1_ptheta,"same");

EmbeddedCanvas ctracknum = new EmbeddedCanvas();
ctracknum.initTimer(1000);
ctracknum.divide(2,1);
ctracknum.cd(0);
ctracknum.draw(h1_numtracksu);
ctracknum.draw(h1_numtracks,"same");
ctracknum.cd(1);
ctracknum.draw(h1_numhitsu);
ctracknum.draw(h1_numhits,"same");

EmbeddedCanvas c_phi = new EmbeddedCanvas();
c_phi.initTimer(1000);
c_phi.divide(2,1);
c_phi.cd(0);
c_phi.draw(h2_phie_vs_phip);
c_phi.cd(1);
c_phi.draw(h1_phidiff);

EmbeddedCanvas c_vz = new EmbeddedCanvas();
c_vz.initTimer(1000);
c_vz.divide(2,1);
c_vz.cd(0);
c_vz.draw(h2_vze_vs_vzp);
c_vz.cd(1);
c_vz.draw(h1_vzdiff);

tabbedPane.add("Electron Kinematis", c_ekin);
tabbedPane.add("Proton Kinematics", c_p1d);
tabbedPane.add("Track info", ctracknum);
tabbedPane.add("Phi", c_phi);
tabbedPane.add("Vz", c_vz);

frame.add(tabbedPane);
frame.setLocationRelativeTo(null);
frame.setVisible(true);

// Define beam lorentz vector and target lorentz vector
LorentzVector target = new LorentzVector(0.0,0.0,0.0,p_mass);

int lineNo = 0;
int run = 0;

// RGF files to open
new File('.', args[0]).eachLine { line ->
    HipoReader reader=new HipoReader();
    reader.open(line);
    
    Event event = new Event();
       
    //new bank definition
    Bank parts = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
    Bank calos = new Bank(reader.getSchemaFactory().getSchema("REC::Calorimeter"));
    Bank chers = new Bank(reader.getSchemaFactory().getSchema("REC::Cherenkov"));
    
    Bank rtpc_tracks = new Bank(reader.getSchemaFactory().getSchema("RTPC::tracks"));
    Bank rtpc_hits = new Bank(reader.getSchemaFactory().getSchema("RTPC::hits"));
    Bank run_config = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));
    
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
        event.read(run_config);
        
        run = run_config.getInt("run", 0); 
        
        if(run>11619 && run<=11656)      beamEnergy=2.14418;
        else if(run>11656)               beamEnergy=10.3894;
        
        LorentzVector   beam = new LorentzVector(0.0,0.0,beamEnergy,beamEnergy);
        
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
                int num_rtpc_tracks = rtpc_tracks.getRows();
                int num_rtpc_hits = rtpc_hits.getRows();
                
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
                            
                    double nu = beamEnergy - e_prime;
                            
                    double mom = vecE.p();
                    double Q2 = 4.0*beamEnergy*e_prime*Math.pow(Math.sin(theta/2.0),2);
                    double W = Math.sqrt(p_mass*p_mass - Q2 + 2*p_mass*nu);
                    //double Q2 = -vecQ2.mass2();
                    //double xB = Q2/(2.0*p_mass*nu);
                    double xB = Q2/(W*W - p_mass*p_mass + Q2);
                            
                    ectot = ecin+ecout;
                             
                    double sf = energy/mom;
                    
                    double e_vz = parts.getFloat("vz",ipart);
                    
                    theta *= 180/Math.PI;
                    phi *= 180/Math.PI;
                    
                    // find new histo mins and maxes
                    if(W < W_min) W_min = W;
                    if(W > W_max && W != Double.NaN) W_max = W;
                    if(Q2 < Q2_min) Q2_min = Q2;
                    if(Q2 > Q2_max) Q2_max = Q2;
                    if(mom < mom_min) mom_min = mom;
                    if(mom > mom_max) mom_max = mom;
                    
                    
                    // fill uncut e- kinematic histos
                    h1_Wu.fill(W);
                    h1_Q2u.fill(Q2);
                    h1_vzeu.fill(e_vz);
                    h1_xBu.fill(xB);
                    h1_thetau.fill(theta);
                    h1_phiu.fill(phi);
                    h1_emomu.fill(mom);
                        
                    // fill electron kinematic histos
                    if(run < 11656 && e_vz > -15 && e_vz < 15 && Q2 > 0.05 && Q2 < 0.1 && W > 0.85 && W < 1.05){
                                h1_W.fill(W);
                                h1_Q2.fill(Q2);
                                h1_vze.fill(e_vz);
                                h1_xB.fill(xB);
                                h1_theta.fill(theta);
                                h1_phi.fill(phi);
                                h1_emom.fill(mom);
                    }
                    else{
                                h1_W.fill(W);
                                h1_Q2.fill(Q2);
                                h1_vze.fill(e_vz);
                                h1_xB.fill(xB);
                                h1_theta.fill(theta);
                                h1_phi.fill(phi);
                                h1_emom.fill(mom);
                    
                    }
                    
                    if (!rtpc_hits) {
                        System.out.println("No RTPC hits!")
                        continue;
                    } 
            
                    // ************************************************************** 
                    // ************************* begin cuts *************************
                    // **************************************************************
                    // Let's look into the rtpc bank       
                    int pads_per_track = 0; 
                                
                        for(int itr = 0; itr < num_rtpc_tracks; itr++){
                                int tid = 0;
                                int _tid = -991;
                                
                                int cid = 0;
                                int _cid = -991;
                                
                                float R_min = 0;
                                float R_max = 0;
                                float R_ = 0;
                                
                                int trkID   = rtpc_tracks.getInt("trkID",itr);
                                
                                float momx   = rtpc_tracks.getFloat("px",itr);
                                float momy   = rtpc_tracks.getFloat("py",itr);
                                float momz   = rtpc_tracks.getFloat("pz",itr);
                                float pmom = Math.sqrt(momx*momx+momy*momy+momz*momz);
                                
                                float numhits   = rtpc_tracks.getInt("nhits",itr);
                                
                                LorentzVector vecP = new LorentzVector(momx,momy,momz,Math.sqrt(pmom*pmom+p_mass*p_mass));
                                
                                float ptheta =  rtpc_tracks.getFloat("theta",itr);
                                double p_phi = vecP.phi();
                                double e_phi = vecE.phi();
                                e_phi *=  180/Math.PI;
                                p_phi *=  180/Math.PI;
                                  
                                double p_vz = rtpc_tracks.getFloat("vz",itr);
                                double delta_vz = e_vz - p_vz;
                                
                                theta = vecE.theta();
                                
                                float p_proton = Math.sqrt(Q2 + (Q2*Q2)/(4*p_mass*p_mass));
                                float ptheta_pred = Math.atan(1.0/((1+beamEnergy/p_mass)*Math.tan(theta/2.0)));
                                ptheta_pred *= 180/Math.PI;
                                
                                float tshift = 0;
                                // REPLACE WITH A FUNCTION TO GET t_shift
                                for(int k = 0; k < num_rtpc_hits; k++){
                                    
                                    float t_shift = rtpc_hits.getFloat("tdiff",k);
                                    tid = rtpc_hits.getInt("trkID",k);
                                    
                                    if(tid == trkID){ tshift = t_shift; break; }
                                }    
                                
                                // fill uncut proton kinematic histos
                                h1_tshiftu.fill(tshift);
                                h1_numtracksu.fill(num_rtpc_tracks);
                                h1_pmomu.fill(pmom);
                                h1_pthetau.fill(ptheta);
                                h1_vzdiffu.fill(e_vz-p_vz);
                                h1_phidiffu.fill(e_phi-p_phi);
                                h1_numhits.fill(numhits);
                                
                                // Make cuts
                                if(e_vz > -25 && e_vz < 25 
                                && p_vz > -25 && p_vz < 25 
                                && delta_vz > -2.5 && delta_vz < 2.5
                                && numhits > 20 
                                && tshift > -200.0 && tshift < 500.0){
                                //&& Q2 > 0.05 && Q2 < 0.1 && W > 0.85 && W < 1.05){
                                    
                                    h1_numtracks.fill(num_rtpc_tracks);
                                
                                    h2_mom.fill(p_proton,pmom);
                                    h1_momdiff.fill(pmom - p_proton);
                                        
                                    h2_ptheta.fill(ptheta_pred, ptheta);
                                    h1_thdiff.fill(ptheta_pred - ptheta); 
                                
                                    for(int k = 0; k < num_rtpc_hits; k++){
                                        //System.out.println("Track TID: " + trkID + ", hit TID: " + tid);
                                                
                                        if(k == 0){ 
                                            _tid = rtpc_hits.getInt("trkID",k);
                                            _cid = rtpc_hits.getInt("id",k);
                                        }
                                        else{
                                            tid = rtpc_hits.getInt("trkID",k);
                                            
                                            if(tid == _tid){
                                                // make sure the hit is apart of the track
                                                // we're looking at
                                                if(tid == trkID){   
                                                    //System.out.println("tid: " + tid + ", hit: " + k + ", cellid: " + cid + ", t_shift: " + tshift);
                                                     
                                                    cid = rtpc_hits.getInt("id",k);
                                                    
                                                    // see if the hit is apart of the 
                                                    // same cellid as the previous hit
                                                    if(cid != _cid && k == 1){pads_per_track+=2;}
                                                    else if(cid != _cid){pads_per_track++;}
                                                    
                                                    // in centimeters
                                                    float _x = rtpc_hits.getFloat("x",k)/10.0; 
                                                    float _y = rtpc_hits.getFloat("y",k)/10.0;
                                                        
                                                    R_ = Math.sqrt(_x*_x + _y*_y);
                                                    if (R_ < R_min) R_min = R_;
                                                    else if (R_ > R_max) R_max = R_;
                                                    
                                                }
                                                // if it's not part the the track we're interested in
                                                // move on to the next hit
                                                else{
                                                    _tid = tid;
                                                }
                                            }
                                            else {
                                                    if(pads_per_track > 0) h1_numhits.fill(pads_per_track);
                                                    _tid = tid;
                                                    pads_per_track = 0;
                                            }
                                        } 
                                    } // end hits loop
                                    
                                    //if(R_min > 2 && R_min < 10
                                    //&& R_max > 2 && R_max < 10){
                                        h1_tshift.fill(tshift);
                                        h1_pmom.fill(pmom);
                                        h1_ptheta.fill(ptheta);
                                        
                                        h2_vze_vs_vzp.fill(p_vz, e_vz);
                                        h2_phie_vs_phip.fill(p_phi, e_phi);
                
                                        h1_vzdiff.fill(e_vz-p_vz);
                                        h1_phidiff.fill(e_phi-p_phi);
                                    //}
                                }
                                
                        }          
                } // end nphe and energy cuts
            } // end charge cut
        } // end particle loop   
    } // end while(event)
    reader.close();
} // end new line

h1_tshift.normalize(h1_tshift);
h1_pmom.normalize(h1_pmom);
h1_ptheta.normalize(h1_ptheta);
h1_vzdiff.normalize(h1_vzdiff);
h1_phidiff.normalize(h1_phidiff);
                                        
h1_Wu.normalize(h1_Wu.integral());
h1_Q2u.normalize(h1_Q2u.integral());
h1_vzeu.normalize(h1_vzeu.integral());
h1_xBu.normalize(h1_xBu.integral());
h1_thetau.normalize(h1_thetau.integral());
h1_phiu.normalize(h1_phiu.integral());
h1_emomu.normalize(h1_emomu.integral());

h1_W.normalize(h1_W.integral());
h1_Q2.normalize(h1_Q2.integral());
h1_vze.normalize(h1_vze.integral());
h1_xB.normalize(h1_xB.integral());
h1_theta.normalize(h1_theta.integral());
h1_phi.normalize(h1_phi.integral());
h1_emom.normalize(h1_emom.integral());
                    
  //  c_ekin.getPad(2).getAxisX().setRange(mom_min,mom_max);
  //  c_ekin.getPad(3).getAxisX().setRange(W_min,W_max);
  //  c_ekin.getPad(4).getAxisX().setRange(Q2_min,Q2_max);
    
    c_ekin.update();

if(run<=11656){
    
    EmbeddedCanvas c_mom = new EmbeddedCanvas();
    c_mom.divide(2,1);
    c_mom.cd(0);
    c_mom.draw(h2_mom);
    c_mom.cd(1);
    c_mom.draw(h1_momdiff);
    
    EmbeddedCanvas c_ptheta = new EmbeddedCanvas();
    c_ptheta.divide(2,1);
    c_ptheta.cd(0);
    c_ptheta.draw(h2_ptheta);
    c_ptheta.cd(1);
    c_ptheta.draw(h1_thdiff);
    
    tabbedPaneP.add("Elastic P Momentum Analysis", c_mom);
    tabbedPaneP.add("Elastic P Theta Analysis", c_ptheta);
    
    pframe.add(tabbedPaneP);
    pframe.setLocationRelativeTo(null);
    pframe.setVisible(true);
    
    c_mom.save("figs/"+run+"/proton/elastic_mom.png");
    c_ptheta.save("figs/"+run+"/proton/elastic_theta.png");
}

ctracknum.save("figs/"+run+"/proton/track_info.png");
c_phi.save("figs/"+run+"/proton/phi.png");
c_vz.save("figs/"+run+"/proton/vz.png");
c_p1d.save("figs/"+run+"/proton/pkinematics.png");

c_ekin.save("figs/"+run+"/electron/ekinematics.png");

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