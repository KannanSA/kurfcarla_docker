# Start FROM the official pre-built image from the QUIP developers.
# This image contains LAMMPS with the ML-QUIP and GAP packages already enabled.
FROM libatomsquip/quip:public

# Set the working directory inside the container
WORKDIR /app

# The base image has python, but we need to install pandas and ase.
# First, copy the requirements file.
COPY requirements.txt .

# Install the python packages using the pip included in the image
RUN pip install --no-cache-dir -r requirements.txt

# Copy all your project files (the python script and project data folder)
COPY . .

# Set the user to root, which is the default for this container
USER root

# FIX: Create a symbolic link from 'lmp_mpi' (the actual command) to 'lmp'.
# This ensures that scripts calling 'lmp' will work correctly in this container.
RUN ln -s /usr/local/bin/lmp_mpi /usr/local/bin/lmp

# Specify the command to run when the container starts.
CMD ["python", "lammps_runner.py"]
