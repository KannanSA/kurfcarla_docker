# Start FROM the official pre-built image from the QUIP developers.
# This image contains LAMMPS with the ML-QUIP and GAP packages already enabled.
FROM libatomsquip/quip:public

# Create a standard, non-root user to run the application.
# This is a security best practice and solves any "run-as-root" issues with mpirun.
RUN useradd --create-home --shell /bin/bash appuser

# --- THIS IS THE KEY FIX ---
# Create a wrapper script that contains the exact, correct command to launch LAMMPS.
# The python script will call this wrapper, which is more robust than calling mpirun directly.
RUN echo '#!/bin/bash' > /usr/local/bin/run_lammps.sh && \
    echo 'exec mpirun -np 1 lmp_mpi "$@"' >> /usr/local/bin/run_lammps.sh && \
    chmod +x /usr/local/bin/run_lammps.sh

# Set the working directory
WORKDIR /app

# Copy requirements file and install python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your project files
COPY . .

# Change ownership of the app directory to the new user
RUN chown -R appuser:appuser /app

# Switch to the non-root user
USER appuser

# Specify the command to run when the container starts.
CMD ["python", "lammps_runner.py"]
