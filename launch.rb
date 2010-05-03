#!/usr/local/bin/ruby
[100,200, 400,600,800,1000,1200,1440].each do |size|
	[1,4,9,16,25,36,49,64].each do |n|
		if (n <= 8)
			20.times { puts "#{size}  #{n}  " + `mpirun -np #{n} ./a.out #{size}` }
		else
			20.times { puts "#{size}  #{n}  " + 
		`mpirun -hostfile nodes -mca plm_rsh_agent rsh -np #{n} ./a.out #{size}` }
		end
	end
end
