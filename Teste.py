import multiprocessing as mp

import yaml
from individuo import Individuo
import random
from gene import Gene

rosetta_conf = None
cont = 0
    
class teste():
    def __init__(self, D, ss3):
        super().__init__()
        with open("config.yaml", 'r') as stream:
            try:
                config = yaml.load(stream,Loader=yaml.Loader)
                self.T1 = config['T1']
                self.T2 = config['T2']
                self.F_lower = config['F_lower']
                self.F_upper = config['F_upper']
                self.F = config['F']
                self.CR = config['CR']
                self.maxAval = config['maxAval']
                self.proteinName = config['proteinName']
                self.protocolo = config['protocolo']
            except yaml.YAMLError as exc:                
                print(exc)
        self.pop = []
        self.especies = []
        self.D = D
        self.ss = ss3
        self.NP = D*5
        self.C_min = (int)(self.NP ** 0.5)
        self.C_max = (int)(self.NP/10)
        self.c = 2
        self.param = 0
        self.m_nmdf = 0
        self.avaliacoes = 0
        
    def call_initializer(self,parametro):
        global rosetta_conf
        rosetta_conf = parametro
    
    if mp.get_start_method() != "spawn":
        def minha_funcao(self,seed, target, psc):
            F_new, CR_new = self.autoAjuste(target.getF(), target.getCR())
            trial = Individuo(F_new, CR_new)
            trial.setGenes(target.getGenes())
            trial.setPose(target.getPose())
            trial.setSS(target.getSS())
            
            r = random.randint(0, self.D-3)
        
            self.gera_trial2(trial, seed, r, CR_new)
            self.fragment_insert(ind=trial, mode=self.get_mode())
            
            if random.uniform(0,1) <= psc:
                # seleção baseada em contato
                self.contact_map(target, trial)
                if trial.getSum_cm() > target.getSum_cm():
                    return(trial)
                else:
                    return(target)
            else:
            # seleção baseada em estrutura
                if trial.getProb_ss() > target.getProb_ss():
                    return(trial)
                else:
                    return(target)
    
    def otimiza(self, seed, especies, psc, conf, pool):        
        results = pool.starmap(self.minha_funcao,  [(seed, target, psc) for target in especies])
        
        return results,290
        
    def autoAjuste(self, F, CR):
        if random.uniform(0, 1) < self.T1:
            F_new = self.F_lower + random.uniform(0,1)*self.F_upper
        else:
            F_new = F
        if random.uniform(0,1) < self.T2:
            CR_new = random.uniform(0,1)
        else:
            CR_new = CR
        return F_new, CR_new
    
    def gera_trial2(self, trial, seed, r, CR_new):
        
        pose = trial.getPose()
        
        # loop para gerar novo indivíduo
        for j in range(r, self.D):
            if self.ss[j] == self.ss[r]:
                if random.uniform(0,1) <= CR_new:
                    trial.genes[j].setPhi(seed.genes[j].getPhi())
                    trial.genes[j].setPsi(seed.genes[j].getPsi())
                    trial.genes[j].setOmega(seed.genes[j].getOmega())
                    trial.genes[j].setSSType(seed.genes[j].getSSType())

                    pose.set_phi((j+1), seed.genes[j].getPhi())
                    pose.set_psi((j+1), seed.genes[j].getPsi())
                    pose.set_omega((j+1), seed.genes[j].getOmega())
            else:
                break

        trial.setPose(pose)
        
    def get_mode(self):
        # sorteia um modelo de fragmento a ser utilizado (3, 3s, 9, 9s)
        rand = random.randint(0, 3)
        mode = ['3', '3s', '9', '9s'][rand]
        return mode
    
    def fragment_insert(self, ind=None, n=100, temp=2.0, mode='3s'):
        global rosetta_conf
        pose = ind.getPose()
                
        mc = rosetta_conf.get_new_mc(pose, rosetta_conf.score3, temp)
        mover = rosetta_conf.get_mer(mode)
        trial = rosetta_conf.get_new_trial_mover(mover, mc)
        mover = trial        

        mc.set_temperature(temp)
        mover.apply(pose)

        pose.assign(mc.lowest_score_pose())
        
        original = rosetta_conf.score3(pose)
        one_more = original is None
        evals = 0

        for _ in range(n):
            mover.apply(pose)
            evals += 1
            if pose.energies().total_energy() < original:
                if not one_more:
                    break
                else:
                    one_more = False
                    
        pose.assign(mc.lowest_score_pose())
        ind.setPose(pose)
        ind.setFitness(mc.lowest_score())
        self.update_angle_from_pose(ind, pose)
        
    def update_angle_from_pose(self, ind, pose):
        genes = []
        rosetta_conf.dssp.apply(pose)    # populates the pose's Pose.secstruct
        ss = pose.secstruct()
        ss_pred = self.ss  # SS predita pelo preditor
        ind.setSS(ss)
        sum_ss = 0
        prob_ss = 0.0
        d = 0
        for k in range(0, self.D):
            gen = Gene()
            gen.setPhi(pose.phi(k + 1))
            gen.setPsi(pose.psi(k + 1))
            gen.setOmega(pose.omega(k + 1))
            gen.setSSType(ss[k])
            genes.append(gen)

            if ss[k] == ss_pred[k]:
                sum_ss += 1

            if ss[k] == 'L':
                prob_ss = prob_ss + (float)(rosetta_conf.probs[d])
            elif ss[k] == 'H':
                prob_ss = prob_ss + (float)(rosetta_conf.probs[d+1])
            else:
                prob_ss = prob_ss + (float)(rosetta_conf.probs[d+2])
            d = d+3
        ind.setGenes(genes)
        ind.setSum_ss(sum_ss)
        ind.setProb_ss(prob_ss)
        
    def contact_map(self, target, trial):
        global rosetta_conf

        dist_trial = 0
        dist_target = 0
        d = 0
        poseTrial = trial.getPose()
        poseTarget = target.getPose()
        dTrial = 0.0
        dTarget = 0.0

        while d < len(rosetta_conf.map):
            res1 = (int)(rosetta_conf.map[d])
            res2 = (int)(rosetta_conf.map[d+1])
            prob = (float)(rosetta_conf.map[d+2])
            type_res1 = 'CB'
            type_res2 = 'CB'
            if rosetta_conf.fasta[res1-1] == 'G':
                type_res1 = 'CA'
            if rosetta_conf.fasta[res2-1] == 'G':
                type_res2 = 'CA'

            dTrial = poseTrial.residue(res1).xyz(type_res1).distance(poseTrial.residue(res2).xyz(type_res2))
            if dTrial <= 8.0:
                dist_trial = dist_trial + prob
            else:
                dist_trial = dist_trial + (prob/dTrial)

            dTarget = poseTarget.residue(res1).xyz(type_res1).distance(poseTarget.residue(res2).xyz(type_res2))
            if dTarget <= 8.0:
                dist_target = dist_target + prob
            else:
                dist_target = dist_target + (prob/dTarget)
            
            d = d+3
        target.setSum_cm(dist_target)
        trial.setSum_cm(dist_trial)

        
        