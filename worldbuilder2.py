#    worldbuilder2.py : University Enrollment Network Simulation
#    Copyright (C) 2020 Philip T. Gressman <gresssman@math.upenn.edu> and Jennifer R. Peck <jpeck1@swarthmore.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import probtools
import universal
import random
import math
import ptracker

def get_parameter(optionsdict,parameter,default):
    value = default
    if parameter in optionsdict:
        value = optionsdict[parameter]
    if '_applied' not in optionsdict:
        optionsdict['_applied'] = {}
    optionsdict['_applied'][parameter] = value
    return value

class BinObject(object):
    def __init__(self):
        self.binlist = []
        self.contents = 0
    def add_bin(self,bottom,top,width):
        self.binlist.append({'bottom' : bottom, 'top' : top, 'width' : width, 'contents' : 0})
    def tally(self,item):
        for key in range(len(self.binlist)):
            if item >= self.binlist[key]['bottom'] and item <= self.binlist[key]['top']:
                self.contents += 1
                self.binlist[key]['contents'] += 1
    def output(self):
        for index,item in enumerate(self.binlist):
            print('%i-%i : %5.2f  ' % (item['bottom'],item['top'],100*item['contents']/self.contents),end=' ')
        print('')
    def reset(self):
        self.contents = 0
        for key in self.binlist:
            self.binlist[key]['contents'] = 0
    def makelist(self,size):
        result = []
        for index in range(size):
            result.append(self.get(index/size))
        return result
    def get(self,fraction):
        total = 0
        for item in self.binlist:
            total += item['width']
        whichindex = None
        cumulative = 0
        for index,item in enumerate(self.binlist):
            if cumulative + item['width'] > fraction*total:
                whichindex = index
                break
            cumulative += item['width']
        theta = (cumulative + self.binlist[index]['width'] - fraction*total) / self.binlist[index]['width']
        return int( 0.5+math.exp( theta * math.log(self.binlist[index]['top']) + (1-theta) * math.log(self.binlist[index]['bottom'])))

##### Create the probability distribution for course selection
##### To quickly assign distinct courses to each student, we assume courses are
##### in collections of 5 that all have the same size. Then we just assign
##### clusters to each student and uniformly randomly determine which courses
##### they take in each assigned cluster.

class University(object):
    def __init__(self,optionsdict):
        self.maximum_section_size = get_parameter(optionsdict,'class_size_limit',150)
        self.contact_upscale_factor = get_parameter(optionsdict,'contact_upscale_factor',1.0)
        self.friendship_contacts = get_parameter(optionsdict,'friendship_contacts',4.0) * self.contact_upscale_factor
        self.academic_contacts = get_parameter(optionsdict,'academic_contacts',4.0)  * self.contact_upscale_factor
        self.broad_social_contacts = get_parameter(optionsdict,'broad_social_contacts',2.0) * self.contact_upscale_factor
        self.department_environmental_contacts = get_parameter(optionsdict,'department_environmental_contacts',4.0) * self.contact_upscale_factor
        self.broad_environmental_contacts = get_parameter(optionsdict,'broad_environmental_contacts',4.0) * self.contact_upscale_factor
        self.residential_neighbors = get_parameter(optionsdict,'residential_neighbors',1.0) * self.contact_upscale_factor

        self.online_transition = get_parameter(optionsdict,'online_transition',30)
        self.residential_rate = 1
        self.social_distancing = get_parameter(optionsdict,'social_distancing', True)
        #if self.online_transition is not False:
            #default_of = 0.5 * (1.0 + self.online_transition / self.maximum_section_size)
        #else:
            #default_of = 1.0
        #self.occupancy_factor = get_parameter(optionsdict,'occupancy_factor',default_of)
        self.crowd_reduction_factor = 1.0
        self.activity_reduction_factor  = 1.0
        if self.online_transition is False:
            self.online_transition = 999999
        self.recitation_rules = [50,20,80] # Classes over 50 have recitations of 20, with 80 student max per assistant
        self.selection_engine = {}
        self.department_selector = {}
        self.fastsubsets = None
        self.test = get_parameter(optionsdict,'test',False)
        self.cohort_data = {}
        self.student_data = {}
        self.class_data = {}
        self.department_data = {}
        self.instructor_data = {}
        self.assistant_data = {}
        self.friendship_data = {}
        self.close_contact_data = {}
        self.absent = {}

        self.class_sizes = BinObject()
        self.class_sizes.add_bin(200,800,0.005)
        self.class_sizes.add_bin(100,199,0.015)
        self.class_sizes.add_bin(50,99,0.08)
        self.class_sizes.add_bin(40,49,0.03)
        self.class_sizes.add_bin(30,39,0.05) ##
        self.class_sizes.add_bin(20,29,0.14)
        self.class_sizes.add_bin(10,19,0.41) ##
        self.class_sizes.add_bin(2,9,0.27)

        self.contact_factor = 1.0
        self._one_time_work()
        self.sectionID = 0
        self.contact_process = None
        self.verbose = get_parameter(optionsdict,'verbose',False)
        self.attendance_bins = get_parameter(optionsdict,'attendance_bins',[[0,0.9],[-100,-10]])
        self.attendance_counts = []
        self.classes = 0

    def get_attendance(self):
        dictresult = {}
        for index,bin in enumerate(self.attendance_bins):
            dictresult['Min Attendance ' + str(bin)] = int(100000 * self.attendance_counts[index][0] / self.classes) / 100000
            dictresult['Attendance ' + str(bin)] = int(100000 * self.attendance_counts[index][1] / self.classes) / 100000
        return dictresult


    def class_contact_ratios(self):
        bins = [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,999999]
        sums = {}
        for binitem in bins:
            sums[binitem] = 0
        for itemid in self.class_data:
            students = len(self.class_data[itemid]['students'])
            contacts = 0.5 * students * (students-1)
            for binitem in bins:
                if students <= binitem and self.class_data[itemid]['type'] != 'plenary':
                    sums[binitem] += contacts
        for binitem in bins:
            sums[binitem] = sums[binitem] / sums[bins[-1]]
        print(sums)

    def _test(self):

        print('===== INFO : Distribution of Class Sizes')
        self.verbose = False
        for index in range(100):
            self.generate()
            final = self.class_data

            for courseno in final:
                classsize = len(final[courseno]['students'])
                self.class_sizes.tally(classsize)
            if index % 10 == 9:
                self.class_sizes.output()

        print('===== CLASS CONTACT RATIOS')
        self.class_contact_ratios()
        print('===== INFO : Contact_Process Expectation')
        total = 0
        self.verbose = False
        for index in range(1):
            self.generate()
            result = self.contact_process.expectation([0,1,2,3,4,5,6])
            factor = 2.0 / 7.0 / (universal.students + universal.instructors)
            for key,value in result.items():
                increment = value * factor
                total += value * factor
                print((key + ' ' * 30)[0:30],value * factor)
            print(total/(index+1))
        self.generate()
        self.mean_field()
        newtotal = 0
        factor = 1.0 / (universal.students + universal.instructors)
        for index in range(7):
            for person in range(universal.students+universal.instructors):
                result = self.contact_process.query(person,index)
                for key,value in result.items():
                    newtotal += value
            print(index,newtotal* factor / (index+1))

    def _threshold_curve(self,*misc,R0=3.8,contacts=19,npi_factor=0.25,class_contacts=8.0):
        classroster = {}
        all_results = []
        broad_contacts = 2.0
        for person in self.student_data:
            classroster[person] = []
            for item in self.student_data[person]['classes']:
                if 'type' not in self.class_data[item] or self.class_data[item]['type'] != 'plenary':
                    classroster[person].append(len(self.class_data[item]['students']))
        sizelist = []
        for item in self.class_data:
            if 'type' not in self.class_data[item] or self.class_data[item]['type'] != 'plenary':
                sizelist.append(len(self.class_data[item]['students']))
        too_big = False
        otherinterv = 7/5 * 2.4 / 7 * npi_factor * 2.0 * R0 / contacts * class_contacts/8.0
        broad_rate = npi_factor * 2.0 * R0 / contacts * broad_contacts
        print('gamma',broad_rate,'base_classroom',otherinterv * universal.in_class_base_rate,'base_friend',otherinterv * universal.in_frnd_base_rate_hi)

        for cutoff in range(5,151):
            total = 0
            students = 0
            for pkey in classroster:
                added = False
                for item in classroster[pkey]:
                    if item < cutoff:
                        total += 1
                        if not added:
                            students += 1
                            added = True
            topsum = 0
            botsum = 0
            for item in sizelist:
                if item < cutoff:
                    comparison_rate = otherinterv * (universal.in_class_base_rate + universal.in_frnd_base_rate_hi / math.sqrt(item))
                    weight  = item / (1 - comparison_rate * item)
                    if comparison_rate * item > 1:
                        too_big = True

                    topsum += weight * comparison_rate * item
                    botsum += item
            repro_eff = (topsum + 0*broad_rate * universal.students/(1-broad_rate)) / (botsum + 0*universal.students/(1-broad_rate)) * (total / students)
            # We don't do broad contacts...
            if not too_big:
                print('?????',cutoff,repro_eff)

            outdict = {'cutoff' : cutoff,'students' : students, 'seats_occupied' : total,
            'seats_per_student' : int(total/students*100)/100}
            repro_eff += broad_rate
            print('Hi',broad_rate)
            if not too_big:
                outdict['R_e'] = repro_eff
            all_results.append(outdict)
            if too_big:
                break
        with open('theory.txt','w') as file:
            file.write(str(all_results))
            file.close()
        with open('theory.csv','w') as file:
            file.write(r'cutoff,$\tilde R$' + '\n')
            for item in all_results:
                file.write(str(item['cutoff']) + ',' + str(item['R_e']) + '\n')
            file.close()

    def _one_time_work(self):
        ### https://www.collegedata.com/college/University-of-Pennsylvania
        ### 2-9 students: 27% of classes
        ### 10-19 students: 41% of classes
        ### 20-29 students: 14% of classes
        ### 30-39 students: 5% of classes
        ### 40-49 students: 3% of classes
        ### 50-99 students: 8% of classes
        ### Over 100 students: 2% of classes

        course_size_pdf = self.class_sizes.makelist(int(universal.classes/5))
        course_size_cdf = []
        cumulative = 0
        for index in range(len(course_size_pdf)):
            cumulative += course_size_pdf[index]
            course_size_cdf.append(cumulative)

        # shortdict = {}
        # for item in course_size_pdf:
        #     if item not in shortdict:
        #         shortdict[item] = 0
        #     shortdict[item] += 5
        # print(shortdict)
        # quit()

        _cohort_course_cumulative = {}
        cohort_course_pdfs = {}
        for cohortno in range(universal.class_cohorts):
            _cohort_course_cumulative[cohortno] = 0
            cohort_course_pdfs[cohortno] = []

        for courseno in range(len(course_size_cdf)):
            parts = probtools.smoothly_partition(course_size_cdf[courseno]/cumulative,universal.class_cohorts)
            for index in range(universal.class_cohorts):
                outvalue = parts[index] * universal.class_cohorts
                cohort_course_pdfs[index].append(outvalue - _cohort_course_cumulative[index])
                _cohort_course_cumulative[index] = outvalue

        self.selection_engine = {}
        for index in range(universal.class_cohorts):
            self.selection_engine[index] = probtools.CustomDistribution(cohort_course_pdfs[index])
        self.cohort_course_pdfs = cohort_course_pdfs
        histogram = []
        for deptno in range(universal.departments):
            histogram.append(math.exp(-2.3*deptno/universal.departments))
        self.department_selector = probtools.CustomDistribution(histogram)

        self.fastsubsets = probtools.FastSubsets(5)

    def generate(self):
        if self.verbose:
            print('===== University Generation: Student Scheduler')
        self.assign_students()
        if self.verbose:
            print('===== University Generation: Department Structure')
        self.assign_departments()
        if self.verbose:
            print('===== University Generation: Instructors and Assistants')
        self.staff_classes()
        if self.verbose:
            print('===== University Generation: Sections')
        self.subdivide_into_sections()
        if self.verbose:
            print('===== University Generation: Friendships')
        self.form_friendships()
        if self.verbose:
            print('===== University Generation: Roommate Selection')
        self.generate_close_contacts_linear()
        if self.verbose:
            print('===== University Generation: Baseline Attendance')
        self.take_attendance()
        if self.verbose:
            print('===== University Generation: Spatiotemporal Identification')
        self.spatiotemporal()
        if self.verbose:
            print('===== University Generation: Academic Contacts')
        self.compoundcontact = ptracker.CompoundContact()
        self.register_academic_contacts(daily_contacts=self.academic_contacts)
        if self.test:
            self.compoundcontact._test(14)
        if self.verbose:
            print('===== University Generation: Environmental Contacts')
        self.register_environmental_contacts(daily_contacts=self.department_environmental_contacts)
        if self.test:
            self.compoundcontact._test(14)
        if self.verbose:
            print('===== University Generation: Friendship Contacts')
        self.register_friendship_contacts(daily_contacts=self.friendship_contacts)
        if self.test:
            self.compoundcontact._test(14)
        if self.verbose:
            print('===== University Generation: Broad Contacts')
        self.register_broad_contacts(daily_contacts=self.broad_environmental_contacts,social_contacts=self.broad_social_contacts)
        if self.test:
            self.compoundcontact._test(14)
        if self.verbose:
            print('===== University Generation: Residential Contacts')
        self.register_residential_contacts(residential_neighbors=self.residential_neighbors)
        if self.test:
            self.compoundcontact._test(14)
        self.classes = len(self.class_data)


    def generate_close_contacts_linear(self):
        self.close_contacts = {}
        size_picker = probtools.DiscreteExponential(self.residential_neighbors / 2)
        for cohort in self.cohort_data:
            iso_roommate = 0
            person_searching = None
            all_roommates = []
            for index,person in enumerate(self.cohort_data[cohort]['students']):
                self.close_contacts[person] = []
                backlog = size_picker.draw()
                for offset in range(backlog):
                    previndex = index - offset - 1
                    if previndex >= 0:
                        prevperson = self.cohort_data[cohort]['students'][previndex]
                        self.close_contacts[prevperson].append(person)
                        self.close_contacts[person].append(prevperson)
        cumulative = 0
        total = 0
        for person in self.close_contacts:
            cumulative += len(self.close_contacts[person])
            total += 1
        for person in range(universal.students):
            if len(self.close_contacts[person]) == 0:
                del self.close_contacts[person]
    def update_query_system(self):
        self.compoundcontact.update()
    def query_transmit(self,person):
        result = self.compoundcontact.query_transmit(person)
        return result
    def query_contacts(self,person,daysback):
        result = self.compoundcontact.query_contacts(person,daysback)
        return result
    def register_departure(self,person):
        self.compoundcontact.absent(person)
        if person not in self.absent:
            self.absent[person] = True
        if person not in self.student_data:
            return
        for classno in self.student_data[person]['classes']:
            old_triple = self.class_data[classno]['attendance']
            new_triple = [min(old_triple[0],old_triple[1]-1),old_triple[1]-1,old_triple[2]]
            self.class_data[classno]['attendance'] = new_triple
            #self.class_data[classno]['intensity_factor'] = self.occupancy_factor * (new_triple[1] + 1) / (new_triple[2]+1)
            self._attendance_bin_me(old_triple,new_triple)

    def register_return(self,person):
        self.compoundcontact.present(person)
        if person in self.absent:
            del self.absent[person]
        if person not in self.student_data:
            return
        for classno in self.student_data[person]['classes']:
            old_triple = self.class_data[classno]['attendance']
            new_triple = [old_triple[0],old_triple[1]+1,max(old_triple[2],old_triple[1]+1)]
            self.class_data[classno]['attendance'] = new_triple
            self._attendance_bin_me(old_triple,new_triple)

    def take_attendance(self):
        self.attendance_counts = []
        for index in range(len(self.attendance_bins)):
            self.attendance_counts.append([0,0,0])
        for classkey,classdat in self.class_data.items():
            attendance = len(classdat['students'])
            classdat['attendance'] = [attendance,attendance,attendance] # min,today,max
            self._attendance_bin_me(None,[attendance,attendance,attendance])

    def _attendance_bin_me(self,old_triple,new_triple):
        if new_triple[2] == 0 or (old_triple is not None and old_triple[2] == 0):
            return
        for index in range(3):
            if old_triple is not None:
                old = old_triple[index]
                oldf = old_triple[index]/old_triple[2]
                oldg = old_triple[index] - old_triple[2]
            else:
                old = -1
                oldf = -1
                oldg = 1
            new = new_triple[index]
            newf = new_triple[index] / new_triple[2]
            newg = new_triple[index] - new_triple[2]
            for binno,bin in enumerate(self.attendance_bins):
                old_use = old
                new_use = new
                if bin[1] <= 1.0:
                    old_use = oldf
                    new_use = newf
                if bin[1] <= 0:
                    old_use = oldg
                    new_use = newg
                if bin[0] <= old_use and old_use <= bin[1]:
                    self.attendance_counts[binno][index] -= 1
                if bin[0] <= new_use and new_use <= bin[1]:
                    self.attendance_counts[binno][index] += 1

    def spatiotemporal(self):
        for person in self.student_data:
            self.student_data[person]['schedule'] = {0 : [], 1: [], 2: [], 3 : [], 4 : [], 5 : [], 6 : []}
            for classid in self.student_data[person]['classes']:
                deptid = self.class_data[classid]['department']
                daylist = self.class_data[classid]['days']
                for day in daylist:
                    self.student_data[person]['schedule'][day].append(deptid)
        for person in self.instructor_data:
            self.instructor_data[person]['schedule'] = {0 : [], 1: [], 2: [], 3 : [], 4 : [], 5 : [], 6 : []}
            for classid in self.instructor_data[person]['classes']:
                deptid = self.class_data[classid]['department']
                daylist = self.class_data[classid]['days']
                for day in daylist:
                    self.instructor_data[person]['schedule'][day].append(deptid)
                for day in range(5):
                    self.instructor_data[person]['schedule'][day].append(deptid)

        for person in self.assistant_data:
            self.student_data[person]['schedule'] = {0 : [], 1: [], 2: [], 3 : [], 4 : [], 5 : [], 6 : []}
            for classid in self.assistant_data[person]['assignments']:
                deptid = self.class_data[classid]['department']
                daylist = self.class_data[classid]['days']
                for day in daylist:
                    self.student_data[person]['schedule'][day].append(deptid)
        for department in self.department_data:
            self.department_data[department]['schedule'] = {0 : {}, 1: {}, 2: {}, 3 : {}, 4 : {}, 5 : {}, 6 : {}}
            self.department_data[department]['around_today'] = {0 : [], 1 : [], 2 : [], 3 : [], 4 : [], 5 : [], 6 : []}
            for classid in self.department_data[department]['classes']:
                if len(self.class_data[classid]['students']) < self.online_transition:
                    daylist = self.class_data[classid]['days']
                    for day in daylist:
                        self.department_data[department]['schedule'][day][classid] = True
                    checkfor = ['students','instructors','assistants']
                    for ptype in checkfor:
                        if ptype in self.class_data[classid]:
                            for person in self.class_data[classid][ptype]:
                                for day in daylist:
                                    self.department_data[department]['around_today'][day].append(person)


    def assign_students(self):
        self.student_data = {}
        self.class_data = {}
        self.cohort_data = {}
        for cohort in range(universal.class_cohorts):
            self.cohort_data[cohort] = {'students' : []}
        for index in range(universal.students):
            self.student_data[index] = {}
            courseload = random.randint(4,5)
            cohort = random.randrange(universal.class_cohorts)
            self.cohort_data[cohort]['students'].append(index)
            for subindex in range(courseload):
                chosen = self.selection_engine[cohort].draw()
                if chosen not in self.student_data[index]:
                    self.student_data[index][chosen] = 1
                else:
                    self.student_data[index][chosen] += 1
            personal_result = [0]*courseload
            added = 0
            for key,value in self.student_data[index].items():
                for subnumber in self.fastsubsets.draw(value):
                    courseno = key * 5 + subnumber
                    personal_result[added] = courseno
                    added += 1
                    if courseno not in self.class_data:
                        self.class_data[courseno] = {'students' : [index]}
                    else:
                        self.class_data[courseno]['students'].append(index)
            self.student_data[index] = { 'cohort' : cohort, 'classes' : personal_result }

    def assign_departments(self):
        self.department_data = {}
        total_classes = 0
        for index in range(universal.departments):
            self.department_data[index] = {'classes' : []}
        for key in self.class_data:
            total_classes += 1
            departmentno = self.department_selector.draw()
            self.department_data[departmentno]['classes'].append(key)
            self.class_data[key]['department'] = departmentno
        assigned_instructors = 0
        firstkey = None
        cumulative_classes = 0
        discretionary_assignments = 0
        numdepts = len(self.department_data)
        for key in self.department_data:
            if firstkey is None:
                firstkey = key
            self.department_data[key]['type'] = 'department'
            self.department_data[key]['days'] = [0,1,2,3,4]
            instructor_share = 1+int((universal.instructors-numdepts) * (len(self.department_data[key]['classes'])+cumulative_classes) / total_classes) - discretionary_assignments
            cumulative_classes += len(self.department_data[key]['classes'])
            self.department_data[key]['instructors'] = list(range(universal.students+assigned_instructors,universal.students+instructor_share+assigned_instructors))
            assigned_instructors += instructor_share
            discretionary_assignments += instructor_share - 1
        if assigned_instructors < universal.instructors:
            self.department_data[firstkey]['instructors'] += list(range(universal.students+assigned_instructors,universal.students+universal.instructors))

    def staff_classes(self):
        self.instructor_data = {}
        self.assistant_data = {}
        section_cap = self.maximum_section_size
        recitation_data = self.recitation_rules
        student_classlist = self.student_data
        class_studentlist = self.class_data
        department_dutylist = self.department_data
        instructor_assignments = {}
        for department in department_dutylist:
            for person in department_dutylist[department]['instructors']:
                if person not in instructor_assignments:
                    instructor_assignments[person] = {'department' : [], 'classes' : []}
                instructor_assignments[person]['department'].append(department)
            classes_in_department = department_dutylist[department]['classes']
            instructor_need_histogram = []
            assistant_need_histogram = []
            histogram_keys    = []
            instructor_supply = department_dutylist[department]['instructors']
            instructor_supply_size = len(instructor_supply)
            needs_subdivision = {}
            for classname in classes_in_department:
                class_size = len(class_studentlist[classname]['students'])
                histogram_keys.append(classname)
                instructor_need_histogram.append(int(1+class_size/section_cap))
                if instructor_need_histogram[-1] > 1:
                    class_studentlist[classname]['type'] = 'plenary'
                if class_size > recitation_data[0]:
                    assistants = int(1 + class_size/recitation_data[2])
                else:
                    assistants = 0
                assistant_need_histogram.append(assistants)
            instructor_assigner = probtools.Histogram(instructor_need_histogram)
            on_instructor_no = 0
            while instructor_assigner.total_occupancy() > 0:
                classname = histogram_keys[instructor_assigner.draw(1,False)]
                assigned_instructor = instructor_supply[on_instructor_no % instructor_supply_size]
                on_instructor_no += 1
                if assigned_instructor not in instructor_assignments:
                    instructor_assignments[assigned_instructor] = {'classes' : []}
                instructor_assignments[assigned_instructor]['classes'].append(classname)
                if 'instructors' not in class_studentlist[classname]:
                    class_studentlist[classname]['instructors'] = []
                class_studentlist[classname]['instructors'].append(assigned_instructor)
            for index,classname in enumerate(histogram_keys):
                if 'assistants' not in class_studentlist[classname]:
                    class_studentlist[classname]['assistants'] = []
                for assistantno in range(assistant_need_histogram[index]):
                    person = random.randrange(0,universal.students)
                    while min(student_classlist[person]['classes']) <= classname:
                        person = random.randrange(0,universal.students)
                    class_studentlist[classname]['assistants'].append(person)
                    if 'assistants' not in self.department_data[department]:
                        self.department_data[department]['assistants'] = []
                    self.department_data[department]['assistants'].append(person)
                    if person not in self.assistant_data:
                        self.assistant_data[person] = {'assignments' : [], 'departments' : []}
                    self.assistant_data[person]['assignments'].append(classname)
                    self.assistant_data[person]['departments'].append(department)
        self.instructor_data = instructor_assignments

    def subdivide_into_sections(self):
        self.sectionID = universal.classes
        needs_subdivision = []
        for key in self.class_data:
            if len(self.class_data[key]['instructors']) > 1:
                needs_subdivision.append([key,self.class_data[key]])
                self.class_data[key]['days'] = []
            else:
                self.class_data[key]['type'] = 'class'
                self.class_data[key]['days'] = probtools.list_select(universal.meeting_schedules)
        needs_recitations = []
        for key,data in needs_subdivision:
            sections_needed = len(data['instructors'])
            size_histogram = [int(1+len(data['students'])/sections_needed)] * sections_needed
            section_sorter = probtools.Histogram(size_histogram)
            section_IDs = list(range(self.sectionID,self.sectionID+sections_needed))
            needs_recitations += section_IDs
            data['sections'] = section_IDs
            for thissectionID in section_IDs:
                instructorID = data['instructors'][thissectionID-self.sectionID]
                self.class_data[thissectionID] = {'type' : 'section', 'department' : data['department'], 'plenary' : key, 'students' : [],
                'instructors' : [instructorID], 'assistants' : []}
                self.class_data[thissectionID]['days'] = probtools.list_select(universal.meeting_schedules)
                self.department_data[data['department']]['classes'].append(thissectionID)
                if key in self.instructor_data[instructorID]['classes']:
                    self.instructor_data[instructorID]['classes'].remove(key)
                if 'plenary' not in self.instructor_data[instructorID]:
                    self.instructor_data[instructorID]['plenary'] = []
                if key not in self.instructor_data[instructorID]['plenary']:
                    self.instructor_data[instructorID]['plenary'].append(key)
                self.instructor_data[instructorID]['classes'].append(thissectionID)

            section_assignment = section_IDs[0]
            for person in data['assistants']:
                self.class_data[section_assignment]['assistants'].append(person)
                if key in self.assistant_data[person]['assignments']:
                    self.assistant_data[person]['assignments'].remove(key)
                self.assistant_data[person]['assignments'].append(section_assignment)
                if 'plenary' not in self.assistant_data[person]:
                    self.assistant_data[person]['plenary'] = []
                self.assistant_data[person]['plenary'].append(key)
                section_assignment +=1
                if section_assignment > section_IDs[-1]:
                    section_assignment = section_IDs[0]

            self.sectionID += sections_needed
            for index,person in enumerate(data['students']):
                newsection = section_IDs[section_sorter.draw(1,False)]
                self.student_data[person]['classes'].remove(key)
                self.class_data[newsection]['students'].append(person)
                if 'plenary' not in self.student_data[person]:
                    self.student_data[person]['plenary'] = []
                self.student_data[person]['plenary'].append(key)
                self.student_data[person]['classes'].append(newsection)
            for thissectionID in section_IDs:
                recitation_groups = probtools.subdivide(self.class_data[thissectionID]['students'],self.recitation_rules[1])
                self.class_data[thissectionID]['recitations'] = []
                for index,group in enumerate(recitation_groups):
                    self.class_data[self.sectionID] = {'type' : 'recitation', 'department' : self.class_data[thissectionID]['department'],
                    'plenary' : self.class_data[thissectionID]['plenary'], 'section' : thissectionID,
                    'assistants' : [], 'students' : []}
                    random_day = random.randrange(0,5)
                    while random_day in self.class_data[thissectionID]['days']:
                        random_day = random.randrange(0,5)
                    self.class_data[self.sectionID]['days'] = [random_day]
                    if len(self.class_data[thissectionID]['assistants']) > 0:
                        myassistant = self.class_data[thissectionID]['assistants'][index % len(self.class_data[thissectionID]['assistants'])]
                        self.class_data[self.sectionID]['assistants'].append(myassistant)
                        self.assistant_data[myassistant]['assignments'].append(self.sectionID)
                    for person in group:
                        self.student_data[person]['classes'].append(self.sectionID)
                        self.class_data[self.sectionID]['students'].append(person)
                    self.class_data[thissectionID]['recitations'].append(self.sectionID)
                    self.sectionID += 1

        sum = [0,0,0]
        size = [0,0,0]
        pairs = [0,0,0]
        meeting_status = []
        for classid in self.class_data:
            if len(self.class_data[classid]['students']) >= self.online_transition:
                self.class_data[classid]['frequency'] = 0
            if 'type' not in self.class_data[classid] or self.class_data[classid]['type'] != 'plenary':
                csize = len(self.class_data[classid]['students'])
                if csize > 0:
                    sum[0] += 1
                    size[0] += csize
                    pairs[0] += csize * (csize-1)
                    if csize < self.online_transition:
                        sum[1] += 1
                        size[1] += csize
                        pairs[1] += csize * (csize - 1)
                status = 1
                if csize >= self.online_transition:
                    status = 0
                meeting_status.append([csize,status,classid])
                self.class_data[classid]['space_upgrade_factor'] = 1.0
        print('+++++ SD_BEFOR: %6i %8i %10i\n+++++ SD_AFTER: %6i %8i %10i' % (sum[0],size[0],pairs[0],sum[1],size[1],pairs[1]))
        meeting_status.sort(reverse=True)
        unoccupied = 0
        for item in meeting_status:
            if item[1] == 1:
                ratio = 1
                new_roomsize = meeting_status[unoccupied][0]
                old_roomsize = item[0]
                unoccupied += 1
                if new_roomsize > 1.5 * old_roomsize and new_roomsize >= 20:
                    if old_roomsize < 10:
                        old_roomsize = 10
                    ratio = old_roomsize / new_roomsize
                    if ratio > 1 or not self.social_distancing:
                        ratio = 1
                    self.class_data[item[2]]['space_upgrade_factor'] = ratio
                    #print(ratio,old_roomsize,new_roomsize)
                sum[2] += ratio
                size[2] += item[0] * ratio
                pairs[2] += item[0] * (item[0]-1) * ratio
        print('+++++ SD_WEIGH: %6i %8i %10i' % (int(sum[2]),int(size[2]),int(pairs[2])))
        if size[1] > 0:
            self.crowd_reduction_factor = pairs[2]/pairs[1]
            self.activity_reduction_factor = size[1]/size[0]
            print('+++++ Activity Reduction: %6.4f  Crowd Reduction: %6.4f' % (self.activity_reduction_factor,self.crowd_reduction_factor))
        for ptype in [self.student_data,self.instructor_data]:
            for person in ptype:
                ptype[person]['physical_days'] = []
                ptype[person]['virtual_days'] = []
                for classid in ptype[person]['classes']:
                    if 'frequency' in self.class_data[classid] and self.class_data[classid]['frequency'] == 0:
                        key = 'virtual_days'
                    else:
                        key = 'physical_days'
                    for day in self.class_data[classid]['days']:
                        ptype[person][key].append(day)
        if not self.social_distancing:
            self.crowd_reduction_factor = 1


    def form_friendships(self):
        groups_formed = 0
        for classid in self.class_data:
            classobj = self.class_data[classid]
            class_size = len(classobj['students'])
            if class_size > 4 and ('type' not in classobj or classobj['type'] != 'plenary'):
                frequency_factor = 1.0
                if 'frequency' in classobj:
                    frequency_factor = classobj['frequency']
                friendgroups = probtools.ordered_subdivide(classobj['students'],int(math.sqrt(class_size)))
                for group in friendgroups:
                    self.friendship_data[groups_formed] = {'type' : 'friendship',
                    'students' : group, 'days' : classobj['days'],'frequency' : frequency_factor}

                    groups_formed += 1
        for classid,classobj in self.department_data.items():
            class_size = len(classobj['instructors'])
            if class_size > 4:
                friendgroups = probtools.ordered_subdivide(classobj['instructors'],int(math.sqrt(class_size)))
                for group in friendgroups:
                    self.friendship_data[groups_formed] = {'type' : 'friendship',
                    'instructors' : group, 'days' : classobj['days']}
                    groups_formed += 1

    def register_academic_contacts(self,*remainder,daily_contacts):
        rate_adjustment = daily_contacts / 4.0 # Calibrated to produce 4.0 on each active day
        rate_adjustment *= 7/5 # Since they only happen on weekdays
        self.compoundcontact.target += 0.5 * daily_contacts
        for itemid in self.class_data:
            classdata = self.class_data[itemid]
            if 'students' in classdata and len(classdata['students']) > 0 and len(classdata['students']) < self.online_transition:
                daylist = classdata['days']
                for day in daylist:
                    context = self.compoundcontact.new_context(day,'academic')
                    context.rate_factor = universal.in_class_base_rate * rate_adjustment * classdata['space_upgrade_factor']
                    context.social_distance_enabled = self.social_distancing
                    context.add_transmitters(classdata['students'])
                    context.add_receivers(classdata['students'])
                    found_instructor = False
                    if 'instructors' in classdata and len(classdata['instructors']) > 0:
                        context.add_transmitters(classdata['instructors'],10)
                        context.add_receivers(classdata['instructors'],5)
                        found_instructor = True
                    if 'assistants' in classdata and len(classdata['assistants']) > 0:
                        if found_instructor:
                            tx = 4
                            dx = 2
                        else:
                            tx = 10
                            dx = 5
                        context.add_transmitters(classdata['assistants'],tx)
                        context.add_receivers(classdata['assistants'],dx)
        for itemid in self.department_data:
            deptdata = self.department_data[itemid]['instructors']
            if len(deptdata) > 1:
                for day in self.department_data[itemid]['days']:
                    context = self.compoundcontact.new_context(day,'departmental')
                    context.rate_factor = universal.in_dept_base_rate * rate_adjustment
                    context.social_distance_enabled = self.social_distancing
                    context.add_transmitters(deptdata)
                    context.add_receivers(deptdata)

    def register_friendship_contacts(self,*rest,daily_contacts):
        context = self.compoundcontact.new_sparse()
        self.compoundcontact.target += 0.5 * daily_contacts
        #context.social_distance_enabled = False  --- Sparse contacts can't be socially distanced anyway
        rate_adjustment = daily_contacts / 4.0
        for itemid in self.friendship_data:
            frienddata = self.friendship_data[itemid]
            ftype = None
            if 'instructors' in frienddata:
                ftype = 'instructors'
                on_day_factor = universal.in_frnd_base_rate_hi * rate_adjustment
                off_day_factor = 0.0
            elif 'students' in frienddata:
                ftype = 'students'
                on_day_factor = universal.in_frnd_base_rate_hi * rate_adjustment
                off_day_factor = universal.in_frnd_base_rate_low * rate_adjustment
            freq_factor = 1.0
            if 'frequency' in frienddata:
                freq_factor = frienddata['frequency']
                on_day_factor *= freq_factor
                off_day_factor *= freq_factor
            if ftype is not None and len(frienddata[ftype]) > 1 and freq_factor > 0:
                factorlist = []
                daylist = frienddata['days']
                for day in range(7):
                    if day in daylist:
                        factorlist.append(on_day_factor)
                    else:
                        factorlist.append(off_day_factor)
                context.add_product_set(frienddata[ftype],frienddata[ftype],factorlist)

    def register_environmental_contacts(self,*rest,daily_contacts):
        rate_adjustment = daily_contacts / 4.0
        self.compoundcontact.target += 0.5 * daily_contacts
        rate_adjustment *= 7/5
        for deptid in self.department_data:
            if 'around_today' in self.department_data[deptid]:
                for day in self.department_data[deptid]['around_today']:
                    if len(self.department_data[deptid]['around_today'][day]) > 1:
                        context = self.compoundcontact.new_context(day,'environmental')
                        context.social_distance_enabled = self.social_distancing
                        context.traceable = False
                        context.rate_factor = universal.in_dept_broad_base_rate * rate_adjustment * self.crowd_reduction_factor
                        context.add_transmitters(self.department_data[deptid]['around_today'][day])
                        context.add_receivers(self.department_data[deptid]['around_today'][day])

    def register_broad_contacts(self,*rest,daily_contacts,social_contacts):
        allpeople = list(range(universal.students+universal.instructors))
        self.compoundcontact.target += 0.5 * (daily_contacts + social_contacts)
        for day in range(7):
            if day <= 4:
                context = self.compoundcontact.new_context(day,'broad')
                context.social_distance_enabled = self.social_distancing
                context.traceable = False
                active_today = []
                for ptype in [self.student_data,self.instructor_data]:
                    for person in ptype:
                        for outday in ptype[person]['physical_days']:
                            if day == outday:
                                active_today.append(person)
                context.rate_factor = daily_contacts / (2.0 * len(active_today))
                context.rate_factor *= 1 / 2.16 * 1.05 # 2.16 expected activities per weekday for students; halfish that many for instructors
                context.rate_factor *= 7/5 # Since they only happen on weekdays
                context.rate_factor *= self.crowd_reduction_factor
                context.add_transmitters(active_today)
                context.add_receivers(active_today)
                print('+++++ Active on Day',day,':',len(active_today))
            if social_contacts > 0:
                context = self.compoundcontact.new_context(day,'broad social')
                context.social_distance_enabled = False
                context.traceable = True
                context.rate_factor = social_contacts / 2 / (universal.students+universal.instructors - 1)
                if day > 4:
                    context.rate_factor *= 2
                context.rate_factor *= 7/9 # Since they're double on weekends
                context.add_transmitters(allpeople)
                context.add_receivers(allpeople)

    def register_residential_contacts(self,*rest,residential_neighbors):
        context = self.compoundcontact.new_permanent()
        self.compoundcontact.target += self.residential_neighbors
        context.rate = 1
        count = 0
        total = 0
        for person in self.close_contacts:
            if len(self.close_contacts[person]) > 0:
                if person in self.close_contacts[person]:
                    print(person,self.close_contacts[person])
                    quit()
                count += 1
                total += len(self.close_contacts[person])
                context.add_product_set([person],self.close_contacts[person],None)
        #print(count,total,total/count,total/20000)
        #context._test()

class HistoryRecord(object):
    def __init__(self):
        self.records = {}
        self.size = 0
        self.all_records = []
        self.information = None
    def reset(self,regenerate):
        if regenerate:
            self.all_records.append(self.records)
        self.records = {}
        self.size = 0
    def record(self,datadict):
        for item in self.records:
            if item in datadict:
                self.records[item].append(datadict[item])
            else:
                self.records[item].append(0)
        for item in datadict:
            if item not in self.records:
                self.records[item] = [0] * (self.size+1)
                self.records[item][-1] = datadict[item]
        self.size += 1
    def output(self,filename):
        with open(filename,'w') as file:
            file.write(str(self.records))
            file.close()
    def output_all(self,filename):
        with open(filename,'w') as file:
            file.write(str(self.compress()))
            file.close()
    def compress(self):
        result = {'_information' : self.information}
        for index,dictionary in enumerate(self.all_records):
            compressed_dict = {}
            for key in dictionary:
                datalist = []
                for datapt in dictionary[key]:
                    if len(datalist) > 0 and datalist[-1][0] == datapt:
                        datalist[-1][1] += 1
                    else:
                        if len(datalist) > 0 and datalist[-1][1] == 1:
                            datalist[-1] = datalist[-1][0]
                        datalist.append([datapt,1])
                if len(datalist) > 0 and type(datalist[-1]) == list and datalist[-1][1] == 1:
                    datalist[-1] = datalist[-1][0]
                compressed_dict[key] = datalist
            for key in compressed_dict:
                if key in result:
                    result[key].append(compressed_dict[key])
                else:
                    if index > 0:
                        result[key] = [[]] * index
                    else:
                        result[key] = []
                    result[key].append(compressed_dict[key])
            for key in result:
                if key not in compressed_dict and key[0] != '_':
                    result[key].append([])
        return result
