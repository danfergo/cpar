require 'PAPI'

def main
  set = PAPI::EventSet::new  #subPapiEvents
  set.add(PAPI::L1_DCM)
  set.add(PAPI::L2_DCM)
  #logToFile

  begin
    set.start #startCounterIncrement
    menu = menu()
    show_results(set) #showresults
  end while menu

  #unsubEvents
end

def show_results(set)
  res = set.stop
  puts ''
  puts 'PAPI results:'
  puts 'L1 DCM: ' << res[0].to_s
  puts 'L2 DCM: ' << res[1].to_s
end

def menu
  puts ''
  puts '==========================='
  puts '1. Multiplication'
  puts '2. Line Multiplication'
  #puts '3. Parallel Multiplication'
  puts '0. End'
  puts 'Selection?: '
  op = gets.chomp.to_i

  if op == 0
    return false
  elsif op < 0 || op > 3
    puts 'Bad option'
    return true
  end

  puts 'Dimensions: lines cols ? '
  line = gets.chomp.split(' ')
  lines = line[0].to_i
  cols = line[1].to_i

  puts ''
  puts 'Working...'
  puts ''

  multiply_matrix(lines, cols, op)

  true
end

def multiply_matrix(l, c, op)

  matrices = init_matrices(l, c)

  start_time = Time.new

  case op
    when 1
      matrices = mult_matrices(matrices, l, c)
    when 2
      matrices = mult_matrices_optimized(matrices, l, c)
    else
      # type code here
  end


  print_results(start_time, op, l, c) #print_results()
  print_result_matrix(matrices[2]) #print_result_matrix()
end

def print_results(start_time, op, m_ar, m_br)
  delta = (Time.new - start_time).round(5)
  performance = 0.0
  #TODO add performance

  puts 'Duration: ' << delta.to_s
end

def mult_matrices(matrices, m_ar, m_br)
  pha = matrices[0]
  phb = matrices[1]
  phc = matrices[2]

  (0..m_ar-1).each do |i|
    (0..m_ar-1).each do |j|
      temp = 0
      (0..m_br-1).each do |k|
        temp += pha[i * m_ar + k] * phb[k * m_br + j]
      end
      phc[i * m_ar + j] = temp
    end
  end
  return pha, phb, phc
  #tentar return matrices
end

def mult_matrices_optimized(matrices, m_ar, m_br)
  pha = matrices[0]
  phb = matrices[1]
  phc = matrices[2]

  (0..m_ar-1).each do |i|
    (0..m_br-1).each do |k|
      (0..m_ar-1).each do |j|
        phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j]
      end
    end
  end
  return pha, phb, phc
  #tentar return matrices
end

def init_matrices(m_ar, m_br)
  pha = []
  phb = []
  phc = []

  (0..m_ar-1).each do |i|
    (0..m_br-1).each do |j|
      pha[i * m_br + j] = 1.0
    end
  end

  (0..m_br-1).each do |i|
    (0..m_ar-1).each do |j|
      phb[i * m_ar + j] = i + 1.0
    end
  end

  (0..m_ar-1).each do |i|
    (0..m_ar-1).each do |j|
      phc[i * m_ar + j] = 0.0
    end
  end

  return pha, phb, phc
end

def print_result_matrix(matrix)
  puts 'Result matrix:'
  i = 0
  matrix.each do |element|
    print element
    print ' '
    if i > 10

      break
    end
    i += 1
  end
  puts ''
end

def papi_test
  puts "Found PAPI #{PAPI::VERSION}"
  puts '-----------'
  set = PAPI::EventSet::new
  puts set.possible
  puts '-----------'
  set.add(PAPI::L1_DCM)
  set.add_named('PAPI_L2_DCM')
  puts set.possible
  set.start
  puts vals = set.stop
  set.start
  set.accum(vals)
  puts vals
  puts set.stop
  puts set.read
  puts set.events
  puts set.read_ts
  puts '-----------'
  set = PAPI::EventSet::new
  puts set.possible(false).length
  if PAPI::COMPONENTS.length > 1
    puts '-----------'
    set = PAPI::EventSet::new
    set.assign_component(PAPI::COMPONENTS[1])
    puts set.possible(false)
  end

end

main
#papi_test
