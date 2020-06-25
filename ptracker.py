#    ptracker.py : Contact Simulation Engine for COVID-19 Transmission Simulation
#    Copyright (C) 2020 Philip T. Gressman <gresssman@math.upenn.edu> and Jennifer R. Peck <jpeck1@swarthmore.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, version 3 of the License
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import random
import probtools
import universal

def dictionary_sum(dict1,dict2):
    if len(dict1) > len(dict2):
        dict3 = dict1
        dict1 = dict2
        dict2 = dict3
    for key in dict1:
        value = dict1[key]
        if key not in dict2:
            dict2[key] = value
        else:
            dict2[key] += value
    return dict2
def dictionary_sum_copy(dict1,dict2):
    result = {}
    for dictN in [dict1,dict2]:
        for key in dictN:
            value = dictN[key]
            if key not in result:
                result[key] = value
            else:
                result[key] += value
    return result

class PersonTracker(object):
    def __init__(self):
        self.ordered_people = {}
        self.person_positions = {}
        self.divider = 0 # Where 'On' Starts; strictly below this is off
        self.divider_memory = 0
        self.total = 0
        self.queue = {}
        self.active = True
        self.switch_to = 0
        self.addition_tasks = {}
    def total_length(self):
        return self.total
    def active_length(self):
        return self.total - self.divider
    def weight(self,person):
        if person not in self.person_positions:
            return 0
        return len(self.person_positions[person])
    def activate(self):
        if self.active is False:
            self.restore()
            self.active = True
            for person in self.queue:
                self.set_state(person,self.queue[person])
            self.queue = {}
    def deactivate(self,switch_to = 0):
        if self.active is True:
            self.save()
            self.active = False
            self.switch_to = switch_to
    def random(self):
        # Returns a random person in the on state, weighted by their multiplicity in the list
        if self.divider == self.total:
            return None
        index = random.randrange(self.divider,self.total)
        return self.ordered_people[index]
    def save(self):
        self.divider_memory = self.divider
    def restore(self):
        self.divider = self.divider_memory
    def add(self,personobj,*remainder,multiplicity=1):
        # add always inserts new people in the "on" state
        if type(personobj) == list:
            personlist = personobj
            isdict = False
            if multiplicity == 1:
                for index,person in enumerate(personobj):
                    self.ordered_people[self.total+index] = person
                    if person not in self.person_positions:
                        self.person_positions[person] = {}
                    self.person_positions[person][self.total+index] = True
                self.total += len(personobj)
                return
        elif type(personobj) == dict:
            personlist = personobj
            isdict = True
        else:
            personlist = [personobj]
            isdict = False
        for person in personlist:
            mymult = multiplicity
            if isdict:
                mymult *= personlist[person]
            if person not in self.person_positions:
                self.person_positions[person] = {}
            for slots in range(mymult):
                self.ordered_people[self.total] = person
                self.person_positions[person][self.total] = True
                self.total += 1
    def _move_to(self,person,new_position):
        newpositions = {}
        new_begin = new_position
        new_end = new_begin + len(self.person_positions[person])
        occupy_target = new_position
        for seat_no in self.person_positions[person]:
            if seat_no >= new_begin and seat_no < new_end:
                # Do nothing; it's already in a spot that you want
                newpositions[seat_no] = True
            else:
                while self.ordered_people[occupy_target] == person:
                    occupy_target += 1
                current_occupant = self.ordered_people[occupy_target]
                # At this point we have seat_no which must be vacated and
                # occupy_target which is the first desired seat not already occupied
                del self.person_positions[current_occupant][occupy_target]
                self.person_positions[current_occupant][seat_no] = True
                newpositions[occupy_target] = True
                self.ordered_people[occupy_target] = person
                self.ordered_people[seat_no] = current_occupant
        self.person_positions[person] = newpositions
    def get_state(self,person):
        if person not in self.person_positions:
            return -1
        for position in self.person_positions[person]:
            if position >= self.divider:
                return 1
            else:
                return 0
    def set_state(self,person,state,require_active=True):
        if person not in self.person_positions:
            return False
        if self.active is False and require_active is True:
            self.queue[person] = state
            return
        mystate = self.get_state(person)
        if mystate == state:
            return
        tomove = len(self.person_positions[person])
        if state == 0:
            self._move_to(person,self.divider)
            self.divider += tomove
            return
        new_divider = self.divider - tomove
        self._move_to(person,new_divider)
        self.divider = new_divider
    def touch(self,person):
        if not self.active:
            self.set_state(person,self.switch_to,False)

class EasyTracker(object):
    def __init__(self):
        self.absenteelist = {}
        self.day = 0
        self.absenteelist[self.day] = {}
    def update(self):
        self.day += 1
        self.absenteelist[self.day] = {}
        if self.day-1 in self.absenteelist:
            for person in self.absenteelist[self.day-1]:
                self.absenteelist[self.day][person] = True
        if self.day - 7 in self.absenteelist:
            del self.absenteelist[self.day-7]
    def absent(self,person):
        self.absenteelist[self.day][person] = True
    def present(self,person):
        if person in self.absenteelist[self.day]:
            del self.absenteelist[self.day][person]
    def poll_absent(self,person,daydiff=0):
        if self.day - daydiff not in self.absenteelist:
            return False
        if person in self.absenteelist[self.day-daydiff]:
            return True
        return False

class SparseContact(object):
    def __init__(self):
        self.roster = EasyTracker()
        self.person_data = {}
        self.pair_data = {}
        self.execution_data = {0 : {}}
        self.pairs = 0
        self.day = 0
        self.parent = None
        self.id = None
    def _set_parent(self,parent,id):
        self.parent = parent
        self.id = id
    def add_product_set(self,transmitlist,receivelist,dayweight):
        for actionlist in [transmitlist,receivelist]:
            for person in actionlist:
                if person not in self.person_data:
                    self.person_data[person] = {}
                    self.person_data[person]['events'] = []
                if self.pairs not in self.person_data[person]['events']:
                    self.person_data[person]['events'].append(self.pairs)
                if self.parent is not None:
                    for day in range(7):
                        self.parent._register(person,self.id,day)
        pairlist = []
        for persona in transmitlist:
            for personb in receivelist:
                if persona != personb:
                    pairlist.append((persona,personb))
        self.pair_data[self.pairs] = [pairlist,dayweight]
        self.pairs += 1
    def update(self):
        self.day += 1
        self.roster.update()
        self.execution_data[self.day] = {}
        if self.day-7 in self.execution_data:
            del self.execution_data[self.day-7]
    def absent(self,person):
        self.roster.absent(person)
    def present(self,person):
        self.roster.present(person)
    def _execute_product(self,index,daydiff = 0):
        if self.day-daydiff in self.execution_data and index in self.execution_data[self.day-daydiff]:
            return self.execution_data[self.day-daydiff][index]
        if self.day-daydiff not in self.execution_data:
            self.execution_data[self.day-daydiff] = {}
        self.execution_data[self.day-daydiff][index] = []
        howmany = probtools.draw(self.pair_data[index][1][(self.day-daydiff)%7]*len(self.pair_data[index][0]))
        for count in range(howmany):
            pairno = random.randrange(0,len(self.pair_data[index][0]))
            self.execution_data[self.day-daydiff][index].append(self.pair_data[index][0][pairno])
        return self.execution_data[self.day-daydiff][index]
    def query_transmit(self,person,day=None):
        if day is not None:
            if day > self.day:
                self.update()
            daydiff = day - self.day
        else:
            daydiff = 0

        if self.roster.poll_absent(person,daydiff) is True:
            return {}
        returndict = {}
        for itemid in self.person_data[person]['events']:
            pairs = self._execute_product(itemid,daydiff)
            for pair in pairs:
                if pair[0] == person and self.roster.poll_absent(pair[1],daydiff) is False:
                    returndict = dictionary_sum(returndict,{pair[1]:1})
        return returndict
    def query_receive(self,person,day=None):
        if day is not None:
            if day > self.day:
                self.update()
            daydiff = day - self.day
        else:
            daydiff = 0
        if self.roster.poll_absent(person,daydiff) is True:
            return {}
        returndict = {}
        for itemid in self.person_data[person]['events']:
            pairs = self._execute_product(itemid,daydiff)
            for pair in pairs:
                if pair[1] == person and self.roster.poll_absent(pair[0],daydiff) is False:
                    returndict = dictionary_sum(returndict,{pair[0]:1})
        return returndict
    def query_contacts(self,person,day=None):
        return dictionary_sum(self.query_transmit(person,day),self.query_receive(person,day))




class PermanentContact(SparseContact):
    def __init__(self):
        super().__init__()
        self.rate = 1
    def _execute_product(self,index,daydiff=0):
        result = []
        for pair in self.pair_data[index][0]:
            result += [pair] * self.rate
        return result
    def _test(self):
        return
        #print(self.pair_data)



class SimpleContact(object):
    def __init__(self,day=None):
        self.transmitters = PersonTracker()
        self.receivers = PersonTracker()
        self.transmit_events = {}
        self.receive_events = {}
        self.contact_events = {}
        self.rate_factor = 1.0
        self.social_distance_enabled = False
        self.effective_factor = 1.0
        self.parent = None
        self.id = None
        self.day = day
        self.previous_day = None
        self.traceable = True
        self.message = ''
    def _set_parent(self,parent,id):
        self.parent = parent
        self.id = id
    def set_rate(self,ratevalue):
        self.rate_factor = ratevalue
    def add_transmitters(self,persobj,*remainder,multiplicity=1):
        self.transmitters.add(persobj,multiplicity)
        if self.parent is not None:
            self.parent._register(persobj,self.id,self.day)
    def add_receivers(self,persobj,*remainder,multiplicity=1):
        self.receivers.add(persobj,multiplicity)
        if self.parent is not None:
            self.parent._register(persobj,self.id,self.day)
    def compute_factor(self):
        if self.social_distance_enabled:
            self.effective_factor = self.rate_factor * (self.transmitters.active_length() + self.receivers.active_length())/(self.transmitters.total_length() + self.receivers.total_length())
        else:
            self.effective_factor = self.rate_factor
    def present(self,person):
        self.transmitters.set_state(person,1)
        self.receivers.set_state(person,1)
    def absent(self,person):
        self.transmitters.set_state(person,0)
        self.receivers.set_state(person,0)
    def initialize(self):
        self.transmitters.deactivate()
        self.receivers.deactivate()
        self.compute_factor()
    def update(self):
        self.transmitters.activate()
        self.transmitters.deactivate()
        self.receivers.activate()
        self.receivers.deactivate()
        self.transmit_events = {}
        self.receive_events = {}
        self.contact_events = {}
        self.compute_factor()
    def _grab_from(self,which_list,weight,sourceindividual):
        rate = self.effective_factor * weight * which_list.active_length()
        howmany = probtools.draw(rate)
        result_dict = {}
        for person in range(howmany):
            whoitis = which_list.random()
            if whoitis != sourceindividual:
                result_dict = dictionary_sum(result_dict,{whoitis : 1})
        return result_dict
    def query_transmit(self,person,day=None):
        if day is not None and day != self.previous_day:
            self.update() # This works because each one is only run once each week
            self.previous_day = day
        if person not in self.transmit_events:
            self.transmit_events[person] = {}
        if self.transmitters.get_state(person) != 1: # Invalid or inactive
            return self.transmit_events[person]
        weight = self.transmitters.weight(person)
        self.transmitters.touch(person)
        new_information = self._grab_from(self.receivers,weight,person)
        for key in new_information:
            value = new_information[key]
            if key not in self.receive_events:
                self.receive_events[key] = {}
            self.receive_events[key] = dictionary_sum(self.receive_events[key],{person : value})
        self.transmit_events[person] = dictionary_sum(self.transmit_events[person],new_information)
        return self.transmit_events[person]
    def query_receive(self,person,day=None):
        if day is not None and day != self.previous_day:
            self.update()
            self.previous_day = day
        if person not in self.receive_events:
            self.receive_events[person] = {}
        if self.receivers.get_state(person) != 1: # Invalid or inactive
            return self.receive_events[person]
        weight = self.receivers.weight(person)
        self.receivers.touch(person)
        new_information = self._grab_from(self.transmitters,weight,person)
        for key in new_information:
            value = new_information[key]
            if key not in self.transmit_events:
                self.transmit_events[key] = {}
            self.transmit_events[key] = dictionary_sum(self.transmit_events[key],{person : value})
        self.receive_events[person] = dictionary_sum(self.receive_events[person],new_information)
        return self.receive_events[person]
    def query_contacts(self,person,day=None):
        if not self.traceable:
            return {}
        if person not in self.contact_events:
            self.contact_events[person] = dictionary_sum_copy(self.query_transmit(person,day),self.query_receive(person,day))
        return self.contact_events[person]



class CompoundContact(object):
    def __init__(self):
        self.simplecontacts = {}
        self.agents = {}
        self.contacts_by_day = {}
        self.contact_count = 0
        self.day = 0
        self.target = 0
    def update(self):
        self.day += 1
    def _test(self,day_range):
        total = 0
        self.day = 0
        print('----- Transmission Test (Target valid with no social distancing)')
        print('----- Transmission Test Target:',self.target)
        for day in range(day_range):
            for person in self.agents:
                cdict = self.query_transmit(person)
                for key in cdict:
                    total += cdict[key]
            print('Day: %3i  ' % (day),'Average: %8.5f' % (total/(day+1)/len(self.agents)))
            self.update()
        return (total/(day+1)/len(self.agents))

    def new_context(self,day,message=''):
        newcontext = SimpleContact(day)
        newcontext._set_parent(self,self.contact_count)
        newcontext.message = message
        self.simplecontacts[self.contact_count] = newcontext
        self.contact_count += 1
        return newcontext
    def new_sparse(self,message=''):
        newcontext = SparseContact()
        newcontext._set_parent(self,self.contact_count)
        newcontext.message = message
        self.simplecontacts[self.contact_count] = newcontext
        self.contact_count += 1
        return newcontext
    def new_permanent(self,message=''):
        newcontext = PermanentContact()
        newcontext._set_parent(self,self.contact_count)
        newcontext.message = message
        self.simplecontacts[self.contact_count] = newcontext
        self.contact_count += 1
        return newcontext
    def _register(self,persobj,id,day):
        if type(persobj) != list and type(persobj) != dict:
            persobj = [persobj]
        for person in persobj:
            if person not in self.agents:
                self.agents[person] = {}
                self.contacts_by_day[person] = {0 : {}, 1 : {}, 2 : {}, 3: {}, 4 : {}, 5 : {}, 6 : {}}
            if day not in self.contacts_by_day[person]:
                self.contacts_by_day[person][day] = {}
            self.agents[person][id] = True
            self.contacts_by_day[person][day][id] = True
    def present(self,person):
        if person not in self.agents:
            return False
        for id in self.agents[person]:
            self.simplecontacts[id].present(person)
    def absent(self,person):
        if person not in self.agents:
            return False
        for id in self.agents[person]:
            self.simplecontacts[id].absent(person)
    def query_transmit(self,person,offsetday = 0):
        result = {}
        daymod = (self.day + offsetday) % 7
        for contactid in self.contacts_by_day[person][daymod]:
            result = dictionary_sum(result,self.simplecontacts[contactid].query_transmit(person,self.day + offsetday))
        return result
    def query_receive(self,person,offsetday = 0):
        result = {}
        daymod = (self.day + offsetday) % 7
        for contactid in self.contacts_by_day[person][daymod]:
            result = dictionary_sum(result,self.simplecontacts[contactid].query_receive(person,self.day + offsetday))
        return result
    def query_contacts(self,person,offsetday = 0):
        result = {}
        daymod = (self.day + offsetday) % 7
        for contactid in self.contacts_by_day[person][daymod]:
            result = dictionary_sum(result,self.simplecontacts[contactid].query_contacts(person,self.day + offsetday))
        return result
