# Use a lightweight official Python image
FROM python:3.12-slim

# Prevent Python from writing .pyc files and buffering stdout/stderr
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set the working directory inside the container
WORKDIR /app

COPY ./DEplots ./DEplots
COPY setup.py .

# Install your project in editable/development mode
RUN pip install .

# Default co
