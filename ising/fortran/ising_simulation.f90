module ising_mod
  implicit none
  integer, allocatable :: spin(:,:)
  integer :: n
  double precision :: T, B, total_energy
contains
  ! Initialize the lattice and compute initial energy
  subroutine initialize(n_in, init_flag, T_in, B_in)
    implicit none
    integer, intent(in) :: n_in
    character(len=*), intent(in) :: init_flag
    double precision, intent(in) :: T_in, B_in
    integer :: i, j
    double precision :: r
    ! Set up lattice size and parameters
    n = n_in
    allocate(spin(n,n))
    T = T_in
    B = B_in
    total_energy = 0.0d0
    ! Set initial spins
    if (init_flag(1:1) == 'r' .or. init_flag(1:1) == 'R') then
       ! Random initial configuration
       do i = 1, n
         do j = 1, n
           call random_number(r)
           if (r < 0.5d0) then
             spin(i,j) =  1
           else
             spin(i,j) = -1
           end if
         end do
       end do
    else if (init_flag(1:1) == 'a' .or. init_flag(1:1) == 'A') then
       ! Aligned initial configuration (all spins +1)
       spin = 1
    else
       print *, "Error: initial configuration flag must be 'random' or 'aligned'"
       stop 1
    end if
    ! Compute initial total energy with periodic boundary conditions
    total_energy = 0.0d0
    ! Sum nearest-neighbor interactions (each pair counted once)
    do i = 1, n
      do j = 1, n-1
        total_energy = total_energy - spin(i,j) * spin(i, j+1)
      end do
    end do
    ! Horizontal wrap-around bonds (last column with first column)
    do i = 1, n
      total_energy = total_energy - spin(i, n) * spin(i, 1)
    end do
    do i = 1, n-1
      do j = 1, n
        total_energy = total_energy - spin(i,j) * spin(i+1, j)
      end do
    end do
    ! Vertical wrap-around bonds (bottom row with top row)
    do j = 1, n
      total_energy = total_energy - spin(n, j) * spin(1, j)
    end do
    ! External field contribution
    do i = 1, n
      do j = 1, n
        total_energy = total_energy - B * spin(i, j)
      end do
    end do
  end subroutine initialize

  ! Perform one Monte Carlo sweep (n*n flip attempts)
  subroutine metropolis_sweep()
    implicit none
    integer :: a, i, j, up, down, left, right, neighbor_sum
    double precision :: r, dE
    do a = 1, n*n
      ! Pick a random site (i,j)
      call random_number(r)
      i = int(r * n) + 1
      call random_number(r)
      j = int(r * n) + 1
      ! Identify neighbors with periodic boundaries
      up    = i - 1; if (up < 1) up = n
      down  = i + 1; if (down > n) down = 1
      left  = j - 1; if (left < 1) left = n
      right = j + 1; if (right > n) right = 1
      neighbor_sum = spin(up,j) + spin(down,j) + spin(i,left) + spin(i,right)
      ! Energy change if this spin flips
      dE = 2.0d0 * spin(i, j) * (neighbor_sum + B)
      if (dE <= 0.0d0) then
        ! Accept flip that lowers energy
        spin(i, j) = - spin(i, j)
        total_energy = total_energy + dE
      else
        ! Accept flip that raises energy with probability exp(-dE/T)
        call random_number(r)
        if (r < exp(-dE / T)) then
          spin(i, j) = - spin(i, j)
          total_energy = total_energy + dE
        end if
      end if
    end do
  end subroutine metropolis_sweep

  ! Save current spin configuration to a text file
  subroutine save_snapshot(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: i, j, unit_id
    open(newunit = unit_id, file = filename, status = 'replace', action = 'write')
    do i = 1, n
      write(unit_id, *) (spin(i, j), j = 1, n)
    end do
    close(unit_id)
  end subroutine save_snapshot
end module ising_mod

program ising_simulation
  use ising_mod
  implicit none
  integer :: arg_count, steps, num_samples, sample_idx, snapshot_idx
  integer :: energy_unit, n_val, nseed, clock, snapshot_count
  integer, allocatable :: sample_steps(:), snapshot_steps(:), seed_array(:)
  double precision :: temp, field
  character(len=100) :: arg_str
  character(len=20)  :: init_flag
  integer :: i

  ! Read command-line arguments
  arg_count = command_argument_count()
  if (arg_count < 6) then
    print *, "Usage: ising_simulation <grid_size> <MC_steps> <num_samples> <temperature> <mag_field> <initial_flag>"
    print *, "Example: ./ising_simulation 32 10000 100 2.5 0.0 random"
    stop
  end if
  call get_command_argument(1, arg_str); read(arg_str, *) n_val
  call get_command_argument(2, arg_str); read(arg_str, *) steps
  call get_command_argument(3, arg_str); read(arg_str, *) num_samples
  call get_command_argument(4, arg_str); read(arg_str, *) temp
  call get_command_argument(5, arg_str); read(arg_str, *) field
  call get_command_argument(6, init_flag)
  if (steps < num_samples) num_samples = steps  ! cap samples if more than steps

  ! Seed RNG using system clock for randomness
  call random_seed(size = nseed)
  if (nseed > 0) then
     allocate(seed_array(nseed))
     call system_clock(count = clock)
     do i = 1, nseed
        seed_array(i) = mod(clock + 37 * i, 1000000000)
     end do
     call random_seed(put = seed_array)
     deallocate(seed_array)
  end if

  ! Initialize lattice and energy
  call initialize(n_val, init_flag, temp, field)

  ! Prepare output for energy samples
  open(newunit = energy_unit, file = 'energy.txt', status = 'replace', action = 'write')
  allocate(sample_steps(num_samples))
  if (num_samples > 1) then
    do i = 0, num_samples-1
       sample_steps(i+1) = 1 + int( dble(i) * dble(steps-1) / dble(num_samples-1) )
    end do
  else
    sample_steps(1) = steps
  end if

  ! Determine which steps to save snapshots (100 frames or fewer if steps < 100)
  snapshot_count = 100
  if (steps < 100) snapshot_count = steps
  allocate(snapshot_steps(snapshot_count))
  if (snapshot_count > 1) then
    do i = 0, snapshot_count-1
       snapshot_steps(i+1) = 1 + int( dble(i) * dble(steps-1) / dble(snapshot_count-1) )
    end do
  else
    snapshot_steps(1) = steps
  end if

  ! Main Monte Carlo loop
  sample_idx   = 1
  snapshot_idx = 1
  do i = 1, steps
     call metropolis_sweep()
     ! Record energy per spin at designated sample points
     if (sample_idx <= num_samples .and. i == sample_steps(sample_idx)) then
        write(energy_unit, *) total_energy / ( dble(n_val) * dble(n_val) )
        sample_idx = sample_idx + 1
     end if
     ! Save snapshot at designated points
     if (snapshot_idx <= size(snapshot_steps) .and. i == snapshot_steps(snapshot_idx)) then
        write(arg_str, '( "snapshot", I3.3, ".txt" )') snapshot_idx-1
        call save_snapshot(trim(arg_str))
        snapshot_idx = snapshot_idx + 1
     end if
  end do

  close(energy_unit)
  deallocate(sample_steps); deallocate(snapshot_steps); deallocate(spin)
end program ising_simulation
