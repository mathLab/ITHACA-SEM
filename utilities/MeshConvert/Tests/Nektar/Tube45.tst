<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nektar prism/tet mesh of pipe intersection with curvature and tet splitting</description>
    <executable>MeshConvert</executable>
    <parameters>-m tetsplit -m jac:list Tube45.rea Tube45.xml:xml:test</parameters>
    <files>
        <file description="Input File">Tube45.rea</file>
        <file description="Surface file">Tube45.hsf</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
