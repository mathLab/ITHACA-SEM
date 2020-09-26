<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Prism Ortho basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s prism -b Ortho_A Ortho_A Ortho_B -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.5 0.0 1.0 0.5 1.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.42088e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">9.23706e-14</value>
        </metric>
    </metrics>
</test>
