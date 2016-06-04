/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob, Romain Prévost

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>
#include <nori/medium.h>
#include <nori/phase.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) {
    m_bvh = new BVH();
}

Scene::~Scene() {
    delete m_bvh;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
	delete m_scene_medium;
}

void Scene::activate() {
    m_bvh->build();

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
        m_sampler->activate();
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

void Scene::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_bvh->addMesh(mesh);
                m_meshes.push_back(mesh);
                if(mesh->isEmitter())
                    m_emitters.push_back(mesh->getEmitter());
            }
            break;
        
        case EEmitter:
			// Add to the background emitter of the scene
			// We know for a fact that there can be only one distant disk in a scene.
			if (static_cast<Emitter*>(obj)->getEmitterType() == EmitterType::EMITTER_DISTANT_DISK || static_cast<Emitter*>(obj)->getEmitterType() == EmitterType::EMITTER_ENVIRONMENT)
				m_bgEmitter = static_cast<Emitter*>(obj);
            m_emitters.push_back(static_cast<Emitter *>(obj));
            break;

        case ESampler:
            if (m_sampler)
                throw NoriException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw NoriException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;
		case EMedium:
			if (m_scene_medium)
				throw NoriException("There can be only one scene medium per scene!");
			m_scene_medium = static_cast<Medium*>(obj);
			break;
        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

std::string Scene::toString() const {
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

    std::string lights;
    for (size_t i=0; i<m_emitters.size(); ++i) {
        lights += std::string("  ") + indent(m_emitters[i]->toString(), 2);
        if (i + 1 < m_emitters.size())
            lights += ",";
        lights += "\n";
    }

    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  meshes = {\n"
        "  %s  }\n"
        "  emitters = {\n"
        "  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(meshes, 2),
        indent(lights,2)
    );
}

Color3f Scene::getBackground(const Ray3f& ray) const
{
	if (m_bgEmitter != nullptr)
	{
		EmitterQueryRecord eRec;
		eRec.ref = ray.o;
		eRec.wi = ray.d;
		return m_bgEmitter->eval(eRec);
	}
	else
		return Color3f(0.0f);
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
