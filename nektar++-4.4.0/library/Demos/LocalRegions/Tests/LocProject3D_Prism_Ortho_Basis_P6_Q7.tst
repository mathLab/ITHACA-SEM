<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Prism Ortho basis P=6 Q=7</description>
    <executable>LocProject3D</executable>
    <parameters>7 1 1 2 6 6 6 7 7 6 0 0 0  1 0 0  1 1 0  0 1 0  0.5 0 1  0.5 1 1</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.42088e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">9.23706e-14</value>
        </metric>
    </metrics>
</test>
