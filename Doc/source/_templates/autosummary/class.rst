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
   {% if item != "__init__"  and item not in inherited_members %}
      ~{{ name }}.{{ item }}
   {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes
   {% for item in attributes %}
   {% if item not in inherited_members %}
   | :data:`{{ item }}`
   {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods_defns %}
   {% if methods %}
   {% for item in methods %}
   {% if item != "__init__"  and item not in inherited_members %}
   .. automethod::  {{ item }}
   {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes_defns %}
   {% if attributes %}
   {% for item in attributes %}
   {% if item not in inherited_members %}
   .. autoattribute::  {{ item }}
   {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}
