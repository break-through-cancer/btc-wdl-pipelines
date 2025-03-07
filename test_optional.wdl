version 1.0

task generate_files {
  input {
    Boolean optional = false
  }
  
  command {
    echo "default file generated" > default.txt
    if ${optional}; then
      echo "optional file generated" > optional.txt
    fi
  }
  
  output {
    File default_output = "default.txt"
    File? optional_output = "optional.txt"
  }
}

workflow simple_pipeline {
  input {
    Boolean optional = true
  }
  
  call generate_files { input: optional = optional }
  
  output {
    File default_output = generate_files.default_output
    File? optional_output = generate_files.optional_output
  }
}
