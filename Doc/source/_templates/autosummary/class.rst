{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
   {% if item != "__init__" %}
      ~{{ name }}.{{ item }}
   {% endif %}
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

   {% block methods_defns %}
   {% if methods %}
   {% for item in methods %}
   {% if item != "__init__" %}
   .. automethod::  {{ item }}
   {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes_defns %}
   {% if attributes %}
   {% for item in attributes %}
   .. autoattribute::  {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
