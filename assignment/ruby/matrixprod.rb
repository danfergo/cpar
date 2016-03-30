require 'PAPI'

class MatrixProduct

  @file
  @set

  def normal
    init

    while menu
    end

    @file.close
  end

  def auto
    init

    #600x600 to 3000x3000 incremented by 400
    mult(600, 3000, 400, 1)

    @file.close
  end

  def mult(min, max, inc, op)
    n = (max-min)/inc
    (0..n).each do |i|
      multiply_matrix(min + i * inc, min + i * inc, op)
    end
  end

  def mult2(min, n, op)
    (0..(n - 1)).each do |i|
      multiply_matrix(min, min, op)
    end
  end

  def init
    @set = PAPI::EventSet::new
    @set.add(PAPI::L1_DCM)
    @set.add(PAPI::L2_DCM)

    #logToFile
    @file = File.open('log.txt', 'w')

    write_headers
  end

  def write_headers
    @file.syswrite('L1 DCM' << ';')
    @file.syswrite('L2 DCM' << ';')
    @file.syswrite('Function ID' << ';')
    @file.syswrite('Rows' << ';')
    @file.syswrite('Cols' << ';')
    @file.syswrite('Duration (seconds)' << ';')
    @file.syswrite('Performance (MFLOPS)' << ';' << "\n")
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

    @set.start
    start_time = Time.new

    case op
      when 1
        matrices = mult_matrices(matrices, l, c)
      when 2
        matrices = mult_matrices_optimized(matrices, l, c)
      else
        # type code here
    end

    papi_res = @set.stop

    print_results(start_time, papi_res, op, l, c) #print_results()
    print_result_matrix(matrices[2]) #print_result_matrix()
  end

  def print_results(start_time, papi_res, op, m_ar, m_br)
    delta = (Time.new - start_time).round(5)
    performance = ((3 * m_ar * m_br * m_br / 1000000) / delta).round(5)

    puts 'Duration: ' << delta.to_s << ' seconds'
    puts 'Prf.: ' << performance.to_s << ' MFLOPS'
    puts 'PAPI results:'
    puts 'L1 DCM: ' << papi_res[0].to_s
    puts 'L2 DCM: ' << papi_res[1].to_s << "\n"

    @file.syswrite(papi_res[0].to_s << ';')
    @file.syswrite(papi_res[1].to_s << ';')
    @file.syswrite(op.to_s << ';')
    @file.syswrite(m_ar.to_s << ';')
    @file.syswrite(m_br.to_s << ';')
    @file.syswrite(delta.to_s << ';')
    @file.syswrite(performance.to_s << ';' << "\n")
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

end

def main
  obj = MatrixProduct.new
  if ARGV[0] == 'auto'
    obj.auto
  else
    obj.normal
  end
end

main
