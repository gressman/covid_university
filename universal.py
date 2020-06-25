#    universal.py : Global Constants for Simulation Scenarios
#    Copyright (C) 2020 Philip T. Gressman <gresssman@math.upenn.edu> and Jennifer R. Peck <jpeck1@swarthmore.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, version 3 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import probtools



students    = 20000
instructors =  2500
classes     =  3750
departments =   120
meeting_schedules = [[1,3],[1,3],[0,2,4],[0,2,4],[0,2]]

class_cohorts = 8

in_class_base_rate = 0.0145
in_dept_base_rate =  0.116
in_frnd_base_rate_hi = 0.18
in_frnd_base_rate_low = in_frnd_base_rate_hi * 0.25
in_dept_broad_base_rate = 0.00158


#_poisson_computer = probtools.Poisson()
#typical_R0 = 4.0
#typical_contacts_per_day = 10.2
#incubation_period = lambda : _poisson_computer.draw(5.2)
#poisson = _poisson_computer.draw
