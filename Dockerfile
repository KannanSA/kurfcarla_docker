# Use the official image with LAMMPS and QUIP pre-installed.
# This image already contains Python, ASE, and a working LAMMPS.
FROM libatomsquip/quip:public

# Create a non-root user for security best practices.
RUN useradd --create-home --shell /bin/bash appuser

# Set the working directory inside the container.
WORKDIR /app

# Copy only the requirements file first to leverage Docker's build cache.
COPY requirements.txt .

# Install the pandas dependency into the existing python environment.
RUN pip install --no-cache-dir -r requirements.txt

# Copy all your project files into the container.
COPY . .

# Give the new user ownership of the files.
RUN chown -R appuser:appuser /app

# Switch to the non-root user.
USER appuser

# Set the default command to run your python script.
CMD ["python", "lammps_runner.py"]