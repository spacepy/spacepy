{{ fullname }}
{{ underline }}

.. automodule:: {{ fullname }}

   {% block functions %}
   {% if functions %}
   .. rubric:: Functions

   .. autosummary::
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: Classes

   .. autosummary::
      :toctree:
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block data %}
   {% if data %}
   .. rubric:: Data

   .. autosummary::
   {% for item in data %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes
   {% for item in attributes %}
   | :data:`{{ item }}`
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: Exceptions

   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions_defns %}
   {% if functions %}
   {% for item in functions %}
   .. autofunction:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block data_defns %}
   {% if data %}
   {% for item in data %}
   .. autodata:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions_defns %}
   {% if exceptions %}

   {% for item in exceptions %}
   .. autoexception:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes_defns %}
   {% if attributes %}
   {% for item in attributes %}
   .. autodata::  {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
