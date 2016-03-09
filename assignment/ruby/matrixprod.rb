

def produtoInterno(v1,v2,col)
	soma = 0.0;
	for i in 0..col
		soma += v1[i] * v2[i];
	end 
end

def main()
	begin
		#------------------MENU----------------
		puts("1. Multiplication");
		puts("2. Line Multiplication");
		puts("Selection?: ")
		op = gets.chomp.to_i
		if(op == 0)
			break;
		end
		puts("Dimensions: lins cols ?")
		line = gets.chomp.split(' ')
		lins = line[0].to_i
		cols = line[1].to_i
		#-------------------------------------- 

		case op
		when 1
		
		when 2

		end

	end while (op != 0)
end
