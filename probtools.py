#    probtools.py : Various probability-related routines for COVID-19 Transmission Simulation
#    Copyright (C) 2020 Philip T. Gressman <gresssman@math.upenn.edu> and Jennifer R. Peck <jpeck1@swarthmore.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import random
import copy
import math

def random_event(p):
    dice = random.random()
    if dice <= p:
        return True
    return False

def random_threshold(outdict):
    dice = random.random()
    result = {}
    for key in outdict:
        value = outdict[key]
        if dice <= value:
            result[key] = value
    return result

def list_permute(mylist):
    listsize = len(mylist)
    selected = dynamicrange(0,listsize,0)
    result = []
    for itemno in range(listsize):
        ordinal = random.randrange(0,listsize-itemno)
        newindex = index_by_order(selected,ordinal,0)
        change_to_state(selected,newindex,1)
        result.append(mylist[newindex])
    return result

def symmetric_subset(k,n): # symmetricly choose k of the elements 0,...,n-1
    selected = dynamicrange(0,n,0)
    result = []
    for itemno in range(k):
        ordinal = random.randrange(0,n-itemno)
        newindex = index_by_order(selected,ordinal,0)
        change_to_state(selected,newindex,1)
        result.append(newindex)
    return result

class FastSubsets(object):
    def __init__(self,pool):
        self.pool = pool
        self.options = {}
        self.selections = {}
        for index in range(0,pool+1):
            self.selections[index] = self._all_symmetric_subsets(index,pool)
            self.options[index] = len(self.selections[index])
    def draw(self,k):
        index = random.randrange(0,self.options[k])
        return self.selections[k][index]
    def _all_symmetric_subsets(self,k,n):
        if k == 0:
            return [[]]
        if k == n:
            return [list(range(0,n))]
        options1 = self._all_symmetric_subsets(k,n-1)
        options2 = self._all_symmetric_subsets(k-1,n-1)
        result = []
        for item in options1:
            result.append(item)
        for item in options2:
            result.append(item + [n-1])
        return result

def dynamicrange(Nl,Nr,state):
    # [Nl,mid,Nr,LeftLn,LeftOcc0,LeftOcc1,RightLn,RightOcc0,RightOcc1,Left,Right]
    if (Nr - Nl) >= 1:
        mid = int((Nr+Nl)/2)
        leftln = mid-Nl
        rightln = Nr-mid
        if state == 0:
            result = [Nl,mid,Nr,leftln,leftln,0,rightln,rightln,0,[],[]]
        else:
            result = [Nl,mid,Nr,leftln,0,leftln,rightln,0,rightln,[],[]]
    return result


def get_state(dynrange,index):
    if index < dynrange[0]:
        return None
    elif index < dynrange[1]:
        if dynrange[3] == dynrange[4]:
            return 0
        elif dynrange[3] == dynrange[5]:
            return 1
        else:
             return get_state(dynrange[9],index)
    elif index < dynrange[2]:
        if dynrange[6] == dynrange[7]:
            return 0
        elif dynrange[6] == dynrange[8]:
            return 1
        else:
            return get_state(dynrange[10],index)
    else:
        return None

def change_to_state(dynrange,index,newstate):
    if index < dynrange[0]:
        return None
    elif index < dynrange[1]:
        side = 0
        offset = 3
        subind = 9
    elif index < dynrange[2]:
        side = 1
        offset = 6
        subind = 10
    else:
        return None
    if dynrange[offset] == dynrange[offset+2-newstate] and dynrange[offset] != 1: # Need to setup substructure because you're breaking totality
        dynrange[subind] = dynamicrange(dynrange[side],dynrange[side+1],1-newstate)
    dynrange[offset+2-newstate] -= 1
    dynrange[offset+1+newstate] += 1
    if dynrange[offset] == dynrange[offset+1+newstate]:
        dynrange[subind] == []
    else:
        change_to_state(dynrange[subind],index,newstate)

# [Nl,mid,Nr,LeftLn,LeftOcc0,LeftOcc1,RightLn,RightOcc0,RightOcc1,Left,Right]
def index_by_order(dynrange,order,state):
    if order < 0:
        return None
    elif order < dynrange[4+state]:
        if dynrange[4+state] == dynrange[3]:
            return dynrange[0] + order
        else:
            return index_by_order(dynrange[9],order,state)
    elif order - dynrange[4+state] < dynrange[7+state]:
        if dynrange[6] == dynrange[7+state]:
            return dynrange[1]+order-dynrange[4+state]
        else:
            return index_by_order(dynrange[10],order-dynrange[4+state],state)
    else:
        return None


def smoothly_partition(p,k):
    if k == 1:
        return [p]
    top = p**k/k
    remains = p - top # Now sums to (k-1)/k
    result = smoothly_partition(remains * k/(k-1),k-1)
    for index in range(len(result)):
        result[index] *= (k-1)/k
    result.append(top)
    return result

def list_select(mylist):
    if len(mylist) == 0:
        return None
    index = random.randrange(0,len(mylist))
    return mylist[index]

class Poisson(object):
    def __init__(self):
        self.CDFlist = []
        self.startI = 8
        self.endI = 1024
        self.data = {}
        self.whole_pool_intensity = [0]
        self.population = []
        self.poolprocs = {}

        thisI = self.startI
        while thisI <= self.endI:
            result = self._createCDF(thisI)
            self.CDFlist.append(result)
            #infostr = '=' * 5 + ' INFO: Poisson %11i -- [%11i - %11i]' % (thisI,result[0],result[0]+len(result[1])-1)
            #print(infostr)
            thisI *= 2

    def _createCDF(self,N):
        startindex = N
        prop_dens = 1.0
        rightlist = []
        atindex = int(N)
        sum = 1.0
        while prop_dens > 1e-20:
            atindex += 1
            prop_dens *= N / atindex
            sum += prop_dens
            rightlist.append(prop_dens)
        atindex = int(N)
        prop_dens = 1.0
        leftlist = []
        while prop_dens > 1e-20 and atindex > 0:
            prop_dens *= atindex/N
            atindex -= 1
            sum += prop_dens
            leftlist.append(prop_dens)
        leftlist.reverse()
        leftlist += [1.0]
        leftlist += rightlist
        factor = 1.0 / sum
        thecdf = []
        runtot = 0.0
        for index in range(len(leftlist)):
            runtot += leftlist[index] * factor
            thecdf.append(runtot)
            leftlist[index] *= factor
        return [atindex,thecdf]
    def _createPDF(self,N,offset=0):  # Generates the PDF for the random variable offset + Poisson(N-offset)
        result1 = self._createCDF(N-offset)
        mypdf = []
        for index in range(result1[0]+offset):
            mypdf.append(0)
        last = 0
        for value in result1[1]:
            mypdf.append(value-last)
            last = value
        return mypdf
    def _draw_from_CDF(self,CDF):
        uniform = random.random()
        lower = -1
        lcuml =  0
        upper = len(CDF) - 1
        rcuml = CDF[-1] # Let's say that rcuml is always greater than *or equal to* uniform
        if uniform > rcuml:
            return len(CDF)
        while upper - lower > 1:
            middle = int((upper + lower)/2)
            mcuml = CDF[middle]
            if mcuml < uniform:
                lower = middle
                lcuml = mcuml
            else:
                upper = middle
                rcuml = mcuml
        return upper
    def draw(self,intensity):
        while intensity > self.endI:
            self.endI *= 2
            result = self._createCDF(self.endI)
            self.CDFlist.append(result)
            #infostr = '=' * 5 + ' INFO: Poisson %11i -- [%11i - %11i]' % (self.endI,result[0],result[0]+len(result[1])-1)
            #print(infostr)
        drawn = 0
        intense_remaining = intensity
        while intense_remaining >= self.startI:
            atI = self.startI
            index = 0
            while intense_remaining > 2*atI:
                atI *= 2
                index += 1
            drawn += self._draw_from_CDF(self.CDFlist[index][1]) + self.CDFlist[index][0]
            intense_remaining -= atI
        repeats = 1 + int(intense_remaining/2)
        ifrac = intense_remaining / repeats
        for rounds in range(repeats):
            incremental_term = math.exp(-ifrac)
            dice = random.random()
            cumulative = incremental_term
            newdraw = 0
            while cumulative < dice:
                newdraw += 1
                incremental_term *= ifrac / newdraw
                cumulative += incremental_term
                if incremental_term < 1e-20:
                    dice = 0
            drawn += newdraw
        return drawn

class AsymmetricProcess(object):
    def __init__(self,itemlist0,itemlist1,intensity_density,poissonproc,parentobj):
        self.parentobj = parentobj
        self.itemlist = [itemlist0,itemlist1]
        self.itemdict = {}
        for index,item in enumerate(itemlist0):
            self.itemdict[item] = (index,0)
        for index,item in enumerate(itemlist1):
            self.itemdict[item] = (index,1)
        self.size = (len(itemlist0),len(itemlist1))
        self.data = { 0 : [dynamicrange(0,self.size[0],0),dynamicrange(0,self.size[1],0)] } # Not queried
        self.queries = { 0 : {}}
        self.intens_dens = intensity_density
        self.poissonproc = poissonproc
        self.state = 0
        self.message = ''
        if type(intensity_density) == list:
            self.current_intensity = intensity_density = intensity_density[0]
        else:
            self.current_intensity = intensity_density
        self.past_states = {}
        self.keep_states = 7
    def reset(self):
        self.data = { 0 : [dynamicrange(0,self.size[0],0),dynamicrange(0,self.size[1],0)] } # Not queried
        self.queries = { 0 : {}}
        self.past_states = {}
    def set_state(self,state):
        if state in self.past_states:
            raise Exception("You can't go back to that state--it's deleted.")
        self.state = state
        if state not in self.data:
            self.data[state] = [dynamicrange(0,self.size[0],0),dynamicrange(0,self.size[1],0)]
            self.queries[state] = {}
        if type(self.intens_dens) == list:
            self.current_intensity = self.intens_dens[state % len(self.intens_dens)]
            if 'intensity_factor' in self.parentobj:
                self.current_intensity *= self.parentobj['intensity_factor']
        self.delete_states()
    def expectation(self):
        value = self.current_intensity * self.size[0] * self.size[1]
        return {self.message : value}
    def delete_states(self):
        todel = {}
        for key in self.queries:
            if key + self.keep_states <= self.state:
                todel[key] = True
                self.past_states[key] = True
        for key in todel:
            del self.queries[key]
            del self.data[key]
    def exhaust(self):
        state = self.state
        avail = {}
        for side in range(2):
            avail[side] = self.data[state][side][4] + self.data[state][side][7]
        howmany = self.poissonproc(self.current_intensity * avail[0] * avail[1])
        for count in range(howmany):
            whichone = random.randrange(0,avail[0])
            whichtwo = random.randrange(0,avail[1])
            pindex = index_by_order(self.data[state][0],whichone,0)
            qindex = index_by_order(self.data[state][1],whichtwo,0)
            pname = self.itemlist[0][pindex]
            qname = self.itemlist[1][qindex]
            pair = [pname,qname]
            for side in range(2):
                if pair[side] not in self.queries[state]:
                    self.queries[state][pair[side]] = {}
                if pair[1-side] not in self.queries[state][pair[side]]:
                    self.queries[state][pair[side]][pair[1-side]] = 0
                self.queries[state][pair[side]][pair[1-side]] += 1
        self.data[state][0] = dynamicrange(0,self.size[0],1) # All open items are settled
        self.data[state][0] = dynamicrange(0,self.size[1],1)
    def query(self,qname,state):
        if state != self.state:
            self.set_state(state)
        qinfo = self.itemdict[qname]
        qindex = qinfo[0]
        qside = qinfo[1]
        if get_state(self.data[state][qside],qindex) == 1:
            if qname not in self.queries[state]:
                return {}
            return self.queries[state][qname]
        change_to_state(self.data[state][qside],qindex,1)
        other_avail = self.data[state][1-qside][4] + self.data[state][1-qside][7]
        if other_avail == 0:
            self.data[state][qside] = dynamicrange(0,self.size[qside],1)
        if qname not in self.queries[state]:
            self.queries[state][qname] = {}
        expectation = self.current_intensity * other_avail
        if expectation < 1:
            self.exhaust()
            if qname not in self.queries[state]:
                return {}
            return self.queries[state][qname]
        howmany = self.poissonproc(expectation)
        for count in range(howmany):
            whichone = random.randrange(0,other_avail)
            index = index_by_order(self.data[state][1-qside],whichone,0)
            pinfo = self.itemlist[1-qside][index]
            pname = pinfo
            if pname not in self.queries[state]:
                self.queries[state][pname] = {}
            if qname not in self.queries[state][pname]:
                self.queries[state][pname][qname] = 0
            if pname not in self.queries[state][qname]:
                self.queries[state][qname][pname] = 0
            self.queries[state][qname][pname] +=1
            self.queries[state][pname][qname] += 1
        return self.queries[state][qname]


class SymmetricProcess(object):
    def __init__(self,itemlist,intensity_density,poissonproc,parentobj):
        self.parentobj = parentobj
        self.itemlist = itemlist
        self.itemdict = {}
        for index,item in enumerate(itemlist):
            self.itemdict[item] = index
        self.size = len(itemlist)
        self.data = { 0 : dynamicrange(0,self.size,0) } # Not queried
        self.queries = { 0 : {}}
        self.intens_dens = intensity_density
        self.poissonproc = poissonproc
        self.state = 0
        if type(intensity_density) == list:
            self.current_intensity = intensity_density = intensity_density[0]
        else:
            self.current_intensity = intensity_density
        self.past_states = {}
        self.keep_states = 7
        self.message = ''
    def reset(self):
        self.data = { 0 : dynamicrange(0,self.size,0) } # Not queried
        self.queries = { 0 : {}}
        self.past_states = {}
    def set_state(self,state):
        if state in self.past_states:
            raise Exception("You can't go back to that state--it's deleted.")
        self.state = state
        if state not in self.data:
            self.data[state] = dynamicrange(0,self.size,0)
            self.queries[state] = {}
        if type(self.intens_dens) == list:
            self.current_intensity = self.intens_dens[state % len(self.intens_dens)]
            if 'intensity_factor' in self.parentobj:
                self.current_intensity *= self.parentobj['intensity_factor']
        self.delete_states()
    def expectation(self):
        value = self.current_intensity * self.size * (self.size - 1)/2
        return { self.message : value}
    def delete_states(self):
        todel = {}
        for key in self.queries:
            if key + self.keep_states <= self.state:
                todel[key] = True
                self.past_states[key] = True
        for key in todel:
            del self.queries[key]
            del self.data[key]
    def exhaust(self):
        state = self.state
        avail = self.data[state][4] + self.data[state][7]
        howmany = self.poissonproc(self.current_intensity * avail * (aval-1)/2)
        for count in range(howmany):
            whichone = random.randrange(0,avail)
            whichtwo = random.randrange(0,avail)
            while whichone == whichtwo:
                whichone = random.randrange(0,avail)
                whichtwo = random.randrange(0,avail)
            pindex = index_by_order(self.data[state],whichone,0)
            qindex = index_by_order(self.data[state],whichtwo,0)
            pname = self.itemlist[pindex]
            qname = self.itemlist[qindex]
            pair = [pname,qname]
            for side in range(2):
                if pair[side] not in self.queries[state]:
                    self.queries[state][pair[side]] = {}
                if pair[1-side] not in self.queries[state][pair[side]]:
                    self.queries[state][pair[side]][pair[1-side]] = 0
                self.queries[state][pair[side]][pair[1-side]] += 1
        self.data[state] = dynamicrange(0,self.size,1) # All open items are settled

    def query(self,qname,state):
        if state != self.state:
            self.set_state(state)
        qindex = self.itemdict[qname]
        if get_state(self.data[state],qindex) == 1:
            if qname not in self.queries[state]:
                return {}
            return self.queries[state][qname]
        change_to_state(self.data[state],qindex,1)
        right_avail = self.data[state][4] + self.data[state][7]
        if qname not in self.queries[state]:
            self.queries[state][qname] = {}
        howmany = self.poissonproc(self.current_intensity * right_avail)
        for count in range(howmany):
            whichone = random.randrange(0,right_avail)
            index = index_by_order(self.data[state],whichone,0)
            pname = self.itemlist[index]
            if pname not in self.queries[state]:
                self.queries[state][pname] = {}
            if qname not in self.queries[state][pname]:
                self.queries[state][pname][qname] = 0
            if pname not in self.queries[state][qname]:
                self.queries[state][qname][pname] = 0
            self.queries[state][qname][pname] +=1
            self.queries[state][pname][qname] += 1
        return self.queries[state][qname]

_global_poisson = Poisson()
draw = _global_poisson.draw

class ContactProcess(object):
    def __init__(self):
        self.actor_processes = {}
        self.processlist = []
        self.poisson_generator = _global_poisson
        self.permanent_contacts = {}
        self.total_permanent_contacts = 0
        self.permanent_contact_rate = 2
    def set_permanent_contacts(self,contactdict):
        self.permanent_contacts = {}
        for person in contactdict:
            thelist = contactdict[person]
            if len(thelist) > 0:
                self.permanent_contacts[person] = {}
                for companion in thelist:
                    self.permanent_contacts[person][companion] = 1
                    self.total_permanent_contacts += 1
        self.total_permanent_contacts = int(self.total_permanent_contacts/2)
    def expectation(self,state):
        totaldict = {'permanent contacts' : self.total_permanent_contacts}
        if type(state) == list:
            for item in state:
                result = self.expectation(item)
                for key in result:
                    value = result[key]
                    if key not in totaldict:
                        totaldict[key] = value
                    else:
                        totaldict[key] += value
            return totaldict
        for process in self.processlist:
            process.set_state(state)
            result =  process.expectation()
            for key in result:
                value = result[key]
                if key not in totaldict:
                    totaldict[key] = value
                else:
                    totaldict[key] += value
        return totaldict
    def actor_length(self):
        return len(self.actor_processes)
    def _register_agents(self,nodelist,whereto):
        for node in nodelist:
            if node not in self.actor_processes:
                self.actor_processes[node] = []
            self.actor_processes[node].append(whereto)
    def register_symmetric(self,nodelist,intensity_densities,message='',parentobj=None):
        if nodelist is None:
            self.register_symmetric(list(self.actor_processes.keys()),intensity_densities,message,parentobj)
            return
        self.processlist.append(SymmetricProcess(nodelist,intensity_densities,self.poisson_generator.draw,parentobj))
        self.processlist[-1].message = message
        whereto = len(self.processlist) - 1
        self._register_agents(nodelist,whereto)
    def register_asymmetric(self,nodelist0,nodelist1,intensity_densities,message='',parentobj=None):
        self.processlist.append(AsymmetricProcess(nodelist0,nodelist1,intensity_densities,self.poisson_generator.draw,parentobj))
        self.processlist[-1].message = message
        whereto = len(self.processlist) - 1
        self._register_agents(nodelist0,whereto)
        self._register_agents(nodelist1,whereto)
    def query(self,person,state):
        result = {}
        if person in self.permanent_contacts:
            for companion in self.permanent_contacts[person]:
                result[companion] = self.permanent_contact_rate
        for process_no in self.actor_processes[person]:
            subresult = self.processlist[process_no].query(person,state)
            for key in subresult:
                value = subresult[key]
                if key not in result:
                    result[key] = value
                else:
                    result[key] += value
        return result


class QuickFind(object):
    def __init__(self,mylist,op,inv):
        self.value = 0
        self.leftQuickFind = None
        self.rightQuickFind = None
        self.terminal = True
        self.op = op
        self.inv = inv
        self.length = len(mylist)
        if self.length > 1:
            midpt = int(self.length/2)
            self.leftQuickFind = QuickFind(mylist[0:midpt],self.op,self.inv)
            self.rightQuickFind = QuickFind(mylist[midpt:self.length],self.op,self.inv)
            self.value = self.op(self.leftQuickFind.value,self.rightQuickFind.value)
            self.terminal = False
        elif self.length == 1:
            self.value = mylist[0]
    def get(self,index):
        if index < 0:
            return None
        elif index >= self.length:
            return None
        if self.terminal:
            return self.value
        if index < self.leftQuickFind.length:
            return self.leftQuickFind.get(index)
        return self.rightQuickFind.get(index-self.leftQuickFind.length)
    def find(self,value=None):
        if value is None:
            value = self.value
        if self.terminal:
            return 0
        if value <= self.leftQuickFind.value:
            return self.leftQuickFind.find(value)
        if value > self.value:
            return None
        result = self.rightQuickFind.find(self.inv(value,self.leftQuickFind.value))
        if result is None:
            return None
        return self.leftQuickFind.length + result
    def change(self,index,value):
        if index < 0 or index >= self.length:
            return None
        if self.terminal:
            self.value = value
            return
        if index < self.leftQuickFind.length:
            self.leftQuickFind.change(index,value)
        else:
            self.rightQuickFind.change(index-self.leftQuickFind.length,value)
        self.value = self.op(self.leftQuickFind.value,self.rightQuickFind.value)

class CustomPDF(QuickFind):
    def __init__(self,pdfobject):
        if type(pdfobject) == dict:
            self.densities = []
            self.labels = []
            for key in pdfobject:
                value = pdfobject[key]
                self.densities.append(value)
                self.labels.append(key)
        elif type(pdfobject) == list:
            self.densities = pdfobject
            self.labels = list(range(len(self.densities)))
        else:
            raise Exception('pdfobject must be a list or a dictionary')
        super().__init__(self.densities,lambda x,y : x+y, lambda x,y : x-y)
    def draw(self):
        dice = random.random() * self.value
        return self.labels[self.find(dice)]

class DiscreteGamma(CustomPDF):
    def __init__(self,mean,cutoff=None):
        distro = {0 : 0}
        increment_factor = (mean-1)/(mean+1)
        basefactor = (1-increment_factor)**2
        stop = False
        atpoint = 1
        cumulative = 0
        mcumulative = 0
        while not stop:
            density = basefactor * atpoint
            #print(density)
            cumulative += density
            mcumulative += density * atpoint
            distro[atpoint] = density
            basefactor *= increment_factor
            atpoint += 1
            if density < 1e-20 or (cutoff is not None and atpoint >= cutoff):
                stop = True
        #print(cumulative,mcumulative)
        super().__init__(distro)
class DiscreteGamma3(CustomPDF):
        def __init__(self,mean,cutoff=None):
            distro = {0 : 0}
            increment_factor = (mean-1)/(mean+2)
            basefactor = (1-increment_factor)**3
            stop = False
            atpoint = 1
            cumulative = 0
            mcumulative = 0
            while not stop:
                density = basefactor * atpoint * (atpoint+1)/2
                #print(density)
                cumulative += density
                mcumulative += density * atpoint
                distro[atpoint] = density
                basefactor *= increment_factor
                atpoint += 1
                if density < 1e-20 or (cutoff is not None and atpoint >= cutoff):
                    stop = True
            #print(cumulative,mcumulative)
            super().__init__(distro)
class DiscreteGammaFull(CustomPDF):
    def __init__(self,mean,shape,cutoff=None):
        distro = {0 : 0}
        increment_factor = (mean-1)/(mean+shape-1)
        basefactor = (1-increment_factor)**shape
        stop = False
        atpoint = 1
        cumulative = 0
        mcumulative = 0
        vcumulative = 0
        while not stop:
            density = basefactor
            for counter in range(shape-1):
                density *= (atpoint+counter) / (counter+1)

            #print(density)
            cumulative += density
            mcumulative += density * atpoint
            vcumulative += density * (atpoint - mean) **2
            distro[atpoint] = density
            basefactor *= increment_factor
            atpoint += 1
            if density < 1e-20 or (cutoff is not None and atpoint >= cutoff):
                stop = True
        #print(cumulative,mcumulative,math.sqrt(vcumulative))
        super().__init__(distro)

class DiscreteExponential(CustomPDF):
    def __init__(self,mean,cutoff=None):
        distro = {0 : 0}
        increment_factor = mean/(mean+1)
        basefactor = (1-increment_factor)
        stop = False
        atpoint = 0
        cumulative = 0
        mcumulative = 0
        while not stop:
            density = basefactor
            #print(density)
            cumulative += density
            mcumulative += density * atpoint
            distro[atpoint] = density
            basefactor *= increment_factor
            atpoint += 1
            if density < 1e-20 or (cutoff is not None and atpoint >= cutoff):
                stop = True
        super().__init__(distro)



class CustomDistribution(object):
    def __init__(self,mylist):
        self.cumulative = {}
        total = 0
        length = 0
        for index,value in enumerate(mylist):
            total += value
            length += 1
            self.cumulative[index] = total
        self.cumulative[-1] = 0
        self.total = total
        self.length = length
        self.decision_tree = {}
        self._populate_decision_tree()

    def _populate_decision_tree(self):
        working_pairs = [(-1,self.length-1)]
        changed = True
        while changed:
            changed = False
            new_pairs = []
            for item in working_pairs:
                if item[1] - item[0] > 1:
                    changed = True
                    midpoint = self._find_midpoint(item)
                    self.decision_tree[item] = midpoint
                    new_pairs += [(item[0],midpoint),(midpoint,item[1])]
            working_pairs = new_pairs
    def _find_midpoint(self,pair):
        bottom = pair[0]
        top = pair[1]
        split = 0.5 * (self.cumulative[bottom] + self.cumulative[top])
        while top - bottom > 1:
            midpt = int((top+bottom)/2)
            midvalue = self.cumulative[midpt]
            if midvalue < split:
                bottom = midpt
            else:
                top = midpt
        if top == pair[0]:
            top += 1
        elif top == pair[1]:
            top -= 1
        return top
    def draw(self):
        dice = random.random() * self.total
        bottom = -1    # We are assuming that the accumulation up to the bottom
        # is strictly less than the dice roll
        top = self.length - 1
        while top - bottom > 1:
            midpt = self.decision_tree[(bottom,top)]
            midvalue = self.cumulative[midpt]
            if midvalue < dice:
                bottom = midpt
            else:
                top = midpt
        return top


class Histogram(object):
    def __init__(self,mylist):
        self.sumobj = QuickFind(mylist,lambda x,y : x+y, lambda x,y : x-y)
        self.maxobj = QuickFind(mylist,lambda x,y : max(x,y), lambda x,y : x)
        self.removal_queue = {}
        self.removal_total = 0
    def total_occupancy(self):
        return self.sumobj.value
    def max_occupancy_bin(self):
        return self.maxobj.find()
    def get_occupancy(self,index):
        return self.sumobj.get(index)
    def set_occupancy(self,index,value):
        self.sumobj.change(index,value)
        self.maxobj.change(index,value)
    def remove(self,index):
        if index not in self.removal_queue:
            self.removal_queue[index] = 0
        self.removal_queue[index] += 1
        self.removal_total + 1
        howmany = self.get_occupancy(index)
        if 4*self.removal_queue[index] > howmany:
            self.set_occupancy(index,max(howmany-self.removal_queue[index],0))
            self.removal_queue[index] = 0

    def draw(self,number=1,replace=True,distro=None):
        if number == 1:
            if distro is None:
                position = random.random() * self.total_occupancy()
                index_drawn = self.sumobj.find(position)
            else:
                chosen = []
                for drawno in range(distro[1]):
                    chosen.append(self.draw(1,True,None))
                chosen.sort()
                index_drawn = chosen[distro[0]]
            if not replace:
                self.remove(index_drawn)
            return index_drawn
        chosen = []
        going = True
        max_attempts = 1
        made_attempts = 0
        to_draw = min(number,self.total_occupancy())
        if to_draw < 1:
            return []
        while going:
            index_drawn = self.draw(1,True,distro)
            if index_drawn in chosen:
                made_attempts += 1
                if made_attempts < max_attempts:
                    chosen = []
                elif made_attempts >= 2 * max_attempts:
                    going = False
            else:
                chosen.append(index_drawn)
                if len(chosen) >= to_draw:
                    going = False
        if not replace:
            for index in chosen:
                self.remove(index)
        return chosen


def subdivide(mylist,target_size):
    listsize = len(mylist)
    if target_size >= listsize:
        return [copy.deepcopy(mylist)]
    histogram = [target_size] * int(1 + len(mylist)/target_size)
    chooser = Histogram(histogram)
    result = {}
    for index,item in enumerate(mylist):
        groupin = chooser.draw(1,False)
        if groupin not in result:
            result[groupin] = []
        result[groupin].append(item)
    list_result = []
    for key in result:
        group = result[key]
        list_result.append(group)
    return list_result

def ordered_subdivide(mylist,target_size):
    listsize = len(mylist)
    if target_size >= listsize:
        return [copy.deepcopy(mylist)]
    howmany = int(1+len(mylist)/target_size)
    list_result = [[]]
    for item in mylist:
        if len(list_result[-1]) >= target_size:
            list_result.append([])
        list_result[-1].append(item)
    return list_result
